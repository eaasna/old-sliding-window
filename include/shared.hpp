#include <seqan3/std/filesystem>
#include <vector>

#include <seqan3/io/sequence_file/input.hpp>

#pragma once

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint8_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct build_arguments
{
    uint32_t window_size{23u};
    uint8_t kmer_size{20u};
    uint8_t threads{1u};
    uint8_t parts{1u};

    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{"./"};
    std::string size{};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool compute_minimiser{false};
    bool compressed{false};
};

struct search_arguments
{
    uint32_t pattern_size{50u};	// this is the size of the sliding window
    // initially using (k,k) minimisers
    uint32_t window_size{17u}; 	// this is the size of the winnowing minimizer
    uint8_t kmer_size{17u};
    uint8_t overlap{49u};	// how much sequential sliding windows i.e patterns overlap
    uint8_t threads{1u};
    uint8_t parts{1u};

    std::filesystem::path query_file{};
    std::filesystem::path ibf_file{};
    std::filesystem::path out_file{"search.out"};
    double tau{0.99};
    double threshold{};
    bool treshold_was_set{false};
    uint8_t errors{0};
    bool compressed{false};
    bool write_time{false};
};
