#include <chrono>
#include <future>
#include <mutex>

// debugging
#include <seqan3/core/debug_stream.hpp>
#include <string>  

// union of result sets
#include <iterator>
#include <set>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <minimiser_model.hpp>
#include <shared.hpp>

class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = default;
    sync_out & operator=(sync_out const &) = default;
    sync_out(sync_out &&) = default;
    sync_out & operator=(sync_out &&) = default;
    ~sync_out() = default;

    sync_out(std::filesystem::path const & path) : file(std::ofstream{path}) {}

    void write(std::string const & data)
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << data;
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};

std::vector<size_t> compute_simple_model(search_arguments const & arguments)
{
    std::vector<size_t> precomp_thresholds;

    if (!arguments.threshold && !do_cerealisation_in(precomp_thresholds, arguments))
    {
// from minimiser_model.cpp
// if k=w then just calculate k-mer lemma threshold
	precomp_thresholds = precompute_threshold(arguments.pattern_size,
                                                  arguments.window_size,
                                                  arguments.kmer_size,
                                                  arguments.errors,
                                                  arguments.tau);

        do_cerealisation_out(precomp_thresholds, arguments);
    }

    return precomp_thresholds;
}

template <typename t>
void load_ibf(t & ibf, search_arguments const & arguments, size_t const part, double & ibf_io_time)
{
    std::filesystem::path ibf_file{arguments.ibf_file};
    ibf_file += "_" + std::to_string(part);
    std::ifstream is{ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    auto start = std::chrono::high_resolution_clock::now();
    iarchive(ibf);
    auto end = std::chrono::high_resolution_clock::now();
    ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <typename t>
inline void do_parallel(t && worker, size_t const num_records, size_t const threads, double & compute_time)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<decltype(std::async(std::launch::async, worker, size_t{}, size_t{}))> tasks;
    size_t const records_per_thread = num_records / threads;

    for (size_t i = 0; i < threads; ++i)
    {
// start/end here show how many records to assign each task
        size_t const start = records_per_thread * i;
        size_t const end = i == (threads-1) ? num_records: records_per_thread * (i+1);
        tasks.emplace_back(std::async(std::launch::async, worker, start, end));
    }

    for (auto && task : tasks)
        task.wait();

    auto end = std::chrono::high_resolution_clock::now();
    compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

template <bool compressed>
void run_program_multiple(search_arguments const & arguments)
{
    constexpr seqan3::data_layout ibf_data_layout = compressed ? seqan3::data_layout::compressed :
                                                                 seqan3::data_layout::uncompressed;
    auto ibf = seqan3::interleaved_bloom_filter<ibf_data_layout>{};

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
                                                std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    size_t const kmer_lemma = arguments.pattern_size + 1u > (arguments.errors + 1u) * arguments.kmer_size ?
                                arguments.pattern_size + 1u - (arguments.errors + 1u) * arguments.kmer_size :
                                0;
    size_t const max_number_of_minimisers = arguments.pattern_size - arguments.window_size + 1;
    std::vector<size_t> const precomp_thresholds = compute_simple_model(arguments);

    auto cereal_worker = [&] ()
    {
        load_ibf(ibf, arguments, 0, ibf_io_time);
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        auto cereal_handle = std::async(std::launch::async, cereal_worker);

        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        std::vector<seqan3::counting_vector<uint16_t>> counts(records.size(),
                                                              seqan3::counting_vector<uint16_t>(ibf.bin_count(), 0));

        auto count_task = [&](size_t const start, size_t const end)
        {
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;

            auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                           seqan3::window_size{arguments.window_size},
                                                           seqan3::seed{adjust_seed(arguments.kmer_size)});

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                (void) id;
                auto & result = counter.bulk_count(seq | hash_view);
                counts[counter_id++] += result;
            }
        };

        cereal_handle.wait();
        do_parallel(count_task, records.size(), arguments.threads, compute_time);

        for (size_t const part : std::views::iota(1u, static_cast<unsigned int>(arguments.parts - 1)))
        {
            load_ibf(ibf, arguments, part, ibf_io_time);
            do_parallel(count_task, records.size(), arguments.threads, compute_time);
        }

        load_ibf(ibf, arguments, arguments.parts - 1, ibf_io_time);
        sync_out synced_out{arguments.out_file};

        auto output_task = [&](size_t const start, size_t const end)
        {
            auto counter = ibf.template counting_agent<uint16_t>();
            size_t counter_id = start;
            std::string result_string{};
            std::vector<uint64_t> minimiser;

            auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                           seqan3::window_size{arguments.window_size},
                                                           seqan3::seed{adjust_seed(arguments.kmer_size)});

            for (auto && [id, seq] : records | seqan3::views::slice(start, end))
            {
                minimiser.clear();
                result_string.clear();
                result_string += id;
                result_string += '\t';

                minimiser = seq | hash_view | seqan3::views::to<std::vector<uint64_t>>;
                counts[counter_id] += counter.bulk_count(minimiser);
                size_t const minimiser_count{minimiser.size()};
                size_t current_bin{0};


                size_t const threshold = arguments.treshold_was_set ?
                                            static_cast<size_t>(minimiser_count * arguments.threshold) :
                                            kmers_per_window == 1 ? kmer_lemma :
                                            precomp_thresholds[std::min(minimiser_count < min_number_of_minimisers ?
                                                                            0 :
                                                                            minimiser_count - min_number_of_minimisers,
                                                                        max_number_of_minimisers -
                                                                            min_number_of_minimisers)] + 2;

                for (auto && count : counts[counter_id++])
                {
                    if (count >= threshold)
                    {
                        result_string += std::to_string(current_bin);
                        result_string += ',';
                    }
                    ++current_bin;
                }
                result_string += '\n';
                synced_out.write(result_string);
            }
        };

        do_parallel(output_task, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
}

// run program single when the index wasn't divided into parts
template <bool compressed>
void run_program_single(search_arguments const & arguments)
{
    constexpr seqan3::data_layout ibf_data_layout = compressed ? seqan3::data_layout::compressed :
                                                                 seqan3::data_layout::uncompressed;
    auto ibf = seqan3::interleaved_bloom_filter<ibf_data_layout>{};

    std::ifstream is{arguments.ibf_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    double ibf_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&] ()
    {
        auto start = std::chrono::high_resolution_clock::now();
        iarchive(ibf);
        auto end = std::chrono::high_resolution_clock::now();
        ibf_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments.out_file};

    // initially kmers_per_window = 1
    size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    // min and max number of minimisers are relevant only for (w,k) minimisers
    size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
                                                std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));

    // used as threshold in case of (k,k)-minimisers when threshold was not set manually
    size_t const kmer_lemma = arguments.pattern_size + 1u > (arguments.errors + 1u) * arguments.kmer_size ?
                                arguments.pattern_size + 1u - (arguments.errors + 1u) * arguments.kmer_size :
                                0;
    size_t const max_number_of_minimisers = arguments.pattern_size - arguments.window_size + 1;
    std::vector<size_t> const precomp_thresholds = compute_simple_model(arguments);

    // [&] captures variables start, end which refer to index of first and last query read in one query file
    auto worker = [&] (size_t const start, size_t const end)
    {
	// sequential calls of bulk_count invalidate the previous result
	auto full_counter = ibf.template counting_agent<uint16_t>();
	auto pre_counter = ibf.template counting_agent<uint16_t>();
	auto post_counter = ibf.template counting_agent<uint16_t>();
        
	std::string result_string{};
	std::set<size_t> result_set; // gather all bins with a hit for each read
        std::vector<uint64_t> minimiser; // minimisers for the whole read
	size_t const threshold = kmer_lemma;

	seqan3::debug_stream << "Pattern size: " << std::to_string(arguments.pattern_size) << '\n';
	seqan3::debug_stream << "Kmer size: " << std::to_string(arguments.kmer_size) << '\n';
	seqan3::debug_stream << "Errors: " << std::to_string(arguments.errors) << '\n';
	seqan3::debug_stream << "Threshold: " << threshold << '\n';

	// all minimisers of a record
        auto hash_view = seqan3::views::minimiser_hash(seqan3::ungapped{arguments.kmer_size},
                                                       seqan3::window_size{arguments.window_size},
                                                       seqan3::seed{adjust_seed(arguments.kmer_size)});


        // && signifies that [id, seq] are r-values?
        // take all records from start to end and loop through them
        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            minimiser.clear();
	    result_set.clear();
            
	    result_string.clear();
            result_string += id;
            result_string += '\t';

            minimiser = seq | hash_view | seqan3::views::to<std::vector<uint64_t>>; // the seq here is the complete read
            
            uint16_t begin; // index of first k-mer in window
            uint16_t end; // index of last k-mer in window 
	    uint16_t pre_begin;
	    uint16_t pre_end;
	    uint8_t len = seq.size(); // length of record

	    //1. find counter result for first pattern minimisers
	    //2. find counter result for minimisers that don't overlap with next pattern
	    //		subtract them from total counter result
	    //3. find counter result for minimisers further from initial pattern
	    //		add them to total counter result


//---------- for first pattern of each read, query the total pattern lenght
            	// TODO: use ranges i.e views instead?
            std::vector<uint64_t> pattern_minimiser(&minimiser[0], 
			    &minimiser[arguments.pattern_size - arguments.kmer_size + 1]);	
		
            auto & result = full_counter.bulk_count(pattern_minimiser); // count occurrences in each bin of all pattern_minimisers
	    std::vector<uint16_t> final_result(result); // changeable vector of initial values
            
	    size_t const minimiser_count{pattern_minimiser.size()}; // not relevant for standard (k,k) minimiser case
            seqan3::debug_stream << "Minimiser count: " << std::to_string(minimiser_count);
	    
	    for (size_t current_bin = 0; current_bin < result.size(); current_bin++)
            {	    
		auto && count = result[current_bin];
                if (count >= threshold)
                {
	        // the result_set is a union of results from all sliding windows of a read
	        result_set.insert(current_bin);

                }
            }
	    
//---------- for sliding window in record do ... for all "middle" patterns
	    for (begin = arguments.pattern_size - arguments.overlap; 
			    begin <= len - arguments.pattern_size; 
			    begin = begin + arguments.pattern_size - arguments.overlap)
            {       
		end = begin + arguments.pattern_size - arguments.kmer_size + 1;

		// store begin and end indices of minimisers in previous pattern
		pre_begin = begin - arguments.pattern_size + arguments.overlap;
		pre_end = end - arguments.pattern_size + arguments.overlap;
		
		// minimisers that were in previous pattern but not in current
		std::vector<uint64_t> pre_pattern_minimiser(&minimiser[pre_begin], &minimiser[begin]);
		auto & pre_result = pre_counter.bulk_count(pre_pattern_minimiser);
		// minimisers that are in current pattern but not in previous
		std::vector<uint64_t> post_pattern_minimiser(&minimiser[pre_end], &minimiser[end]);
		auto & post_result = post_counter.bulk_count(post_pattern_minimiser);
                    
		// TODO: is the new pattern the expected length?
		size_t const new_count{pre_pattern_minimiser.size() + post_pattern_minimiser.size()}; 
		seqan3::debug_stream << "New minimiser count: " << std::to_string(new_count) << '\n';

		for (size_t current_bin = 0; current_bin < result.size(); current_bin++)
		//for (auto && count : result)
                {
		    
		    final_result[current_bin] = result[current_bin] - pre_result[current_bin] + post_result[current_bin];
		    auto && count = final_result[current_bin];
                    if (count >= threshold)
                    {
			// the result_set is a union of results from all sliding windows of a read
			result_set.insert(current_bin);

                    }
                }

            }

//---------- for last pattern in window
// TODO: handle last window separately for ease
/*
	    // last window might have a smaller overlap to make sure that the end of the sequence is covered
	    if (begin + arguments.overlap - 1 < len)
	    {
		// TODO: write a function to do this
	        //returns a reference to the storage of the original vector without copy?
	        std::vector<uint64_t> pattern_minimiser(&minimiser[len - arguments.pattern_size], 
				&minimiser[len - arguments.kmer_size + 1]);
		
		
                seqan3::debug_stream << "Last nr of minimisers: " << pattern_minimiser.size() << 
			'\n' << "Pattern length: " << std::to_string(arguments.pattern_size) << 
			'\n' << "Minimiser length: " << std::to_string(arguments.window_size) << '\n';	
		

		// this counter is the ibf counting agent that gets fed a set of minimisers at a time
                auto & result = counter.bulk_count(pattern_minimiser); // minimiser is a vector of hash values
                size_t const minimiser_count{pattern_minimiser.size()};
                
			
                seqan3::debug_stream << "REsults size: " << std::to_string(result.size()) << '\n'; 
			//<< "" << std::to_string() << << 
	

		size_t current_bin{0};
		for (auto && count : result)
                {
                    if (count >= threshold)
                    {
			// the result_set is a union of results from all sliding windows of a read
			result_set.insert(current_bin);

                    }
                    ++current_bin;
                }
            }
*/
//---------- write out union of bin hits over all patterns
            for (auto bin : result_set)
            {
	        result_string += std::to_string(bin);
	        result_string += ',';
            }
            
	    result_string += '\n';
            synced_out.write(result_string);
        }
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL<<20)*10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::cpp20::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "IBF I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed
                    << std::setprecision(2)
                    << ibf_io_time << '\t'
                    << reads_io_time << '\t'
                    << compute_time;
    }
}

void raptor_search(search_arguments const & arguments)
{
    if (arguments.parts == 1)
    {
        if (arguments.compressed)
            run_program_single<true>(arguments);
        else
            run_program_single<false>(arguments);
    }
    else
    {
        if (arguments.compressed)
            run_program_multiple<true>(arguments);
        else
            run_program_multiple<false>(arguments);
    }

    return;
}
