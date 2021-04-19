#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/std/filesystem>

#include "shared.h"

// ---------------------- Methods for finding the representative k-mer set of the query sequences ----------------------


// take a vector of all k-mers in the sequence and return all k-mers in a window
// window is specified by indices of first and last k-mers in the window
// TODO:
// if the overlap o = w-1 it would be more efficient to use std::deque instead of creating a new vector for each window
// a deque is used in raptor's minimiser_model.cpp but it only allows to pop or push one element at a time
// so would it be more efficent to create a new vector each time?
void find_window(std::vector<uint64_t> all_values, uint64_t begin, uint64_t end)
{
	seqan3::debug_stream << begin << '\t' << end << '\n';
	auto window_values = all_values | 
		seqan3::views::slice(begin, end);
	seqan3::debug_stream << window_values << '\n';
}

// take a vector of all k-mers in the sequence and the sequence length
// return all windows
void find_all_windows(std::vector<uint64_t> all_values, uint8_t len, search_arguments args) 
{
	uint64_t begin; // index of first k-mer in window
	uint64_t end; // index of last k-mer in window
	seqan3::debug_stream << "Text length: " << len << '\n';	
	// get hash values for all windows
	for (begin = 0; begin <= len - args.window; begin = begin + args.window - args.overlap){
		end = begin + args.window - args.kmer_length + 1;
		find_window(all_values, begin, end);
	}
	// last window might have a smaller overlap to make sure that the end of the sequence is covered	
	if (begin + args.overlap - 1 < len){
		find_window(all_values, len - args.window, len - args.window + args.kmer_length);
	}
}

// take a sequence
// return all (k,k)-minimisers of the sequence	
void find_all_kmers(auto text, search_arguments args) 
{
	// The default seed for minimiser_hash
	uint64_t seed = 0x8F3F73B5CF1C9ADE;
	
	// (k, k)-minimizer is used
	std::vector<uint64_t> all_values = text | 
		seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{args.kmer_length}}, 
				seqan3::window_size{args.kmer_length}) | // (k,k)-minimizer
		std::views::transform([seed] (uint64_t i) {return i ^ seed; }) | // use XOR to get back to hash values
		seqan3::views::to<std::vector<uint64_t>>; // re-transform into a distinct container
		
	uint8_t len = text.size(); // text length	
	find_all_windows(all_values, len, args);

	// print all k-mers
	seqan3::debug_stream << all_values << '\n';	
}
