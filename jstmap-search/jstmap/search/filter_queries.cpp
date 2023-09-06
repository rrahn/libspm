// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <omp.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/search/filter_queries.hpp>

namespace jstmap
{

std::tuple<size_t, uint8_t, seqan3::interleaved_bloom_filter<>> load_index(std::filesystem::path const & index_path)
{
    std::ifstream instr{index_path};
    cereal::BinaryInputArchive inarch{instr};
    // Load the bin size used for the jst partitioning.
    size_t bin_size{};
    uint8_t kmer_size{};
    inarch(bin_size);
    inarch(kmer_size);
    // Load the corresponding ibf.
    seqan3::interleaved_bloom_filter<> ibf{};
    ibf.serialize(inarch);

    return std::tuple{bin_size, kmer_size, std::move(ibf)};
}

std::pair<size_t, std::vector<search_queries_type>>
filter_queries(std::vector<search_query> const & queries, search_options const & options)
{
    using bucket_list_t = std::vector<search_queries_type>;
    auto [bin_size, kmer_size, ibf] = load_index(options.index_input_file_path);
    log_debug("IBF bin_size:", bin_size);
    log_debug("IBF kmer_size:", kmer_size);
    log_debug("IBF bin_count:", ibf.bin_count());

    bucket_list_t read_bucket_list{};
    read_bucket_list.resize(ibf.bin_count());

    // size_t const kmers_per_window = arguments.window_size - arguments.kmer_size + 1;
    // size_t const kmers_per_pattern = arguments.pattern_size - arguments.kmer_size + 1;
    // size_t const min_number_of_minimisers = kmers_per_window == 1 ? kmers_per_pattern :
    //                                             std::ceil(kmers_per_pattern / static_cast<double>(kmers_per_window));
    // size_t const kmer_lemma = arguments.pattern_size + 1u > (arguments.errors + 1u) * arguments.kmer_size ?
    //                             arguments.pattern_size + 1u - (arguments.errors + 1u) * arguments.kmer_size :
    //                             0;

    auto counting_agent = ibf.template counting_agent<uint16_t>();
    std::vector<bucket_list_t> thread_local_buffer{};
    thread_local_buffer.resize(options.thread_count, read_bucket_list);

    #pragma omp parallel for num_threads(options.thread_count) shared(thread_local_buffer, queries) firstprivate(counting_agent, kmer_size) schedule(dynamic)
    for (size_t query_idx = 0; query_idx < queries.size(); ++query_idx)
    {
        log_debug("IBF counting query:", query_idx);
        // kmer-lemma:
        size_t const thread_id = omp_get_thread_num();
        search_query query = queries[query_idx];
        size_t const query_size = std::ranges::size(query.value().sequence());
        size_t const error_count = std::floor(query_size * options.error_rate);
        size_t const kmer_threshold = query_size - kmer_size + 1 - (error_count * kmer_size);
        log_debug("IBF kmer_threshold:", kmer_threshold);

        // Counting:
        auto hashed_seq = query.value().sequence() | seqan3::views::kmer_hash(seqan3::ungapped{kmer_size});
        std::vector hashes(std::ranges::begin(hashed_seq), std::ranges::end(hashed_seq));
        auto & bin_counts = counting_agent.bulk_count(hashes);
        // log_debug("IBF:bin_counts", bin_counts);

        // Bin assignment:
        for (size_t bin_idx = 0; bin_idx < bin_counts.size(); ++bin_idx)
            if (bin_counts[bin_idx] >= kmer_threshold) {
                thread_local_buffer[thread_id][bin_idx].push_back(query);
            }
    }

    // Reduce the local buffer and move them to the final bin vector.
    std::ranges::for_each(thread_local_buffer, [&](auto & local_bucket_list) {
        size_t bucket_idx = 0;
        std::ranges::for_each(local_bucket_list, [&] (auto const & source_bucket) {
            search_queries_type & target_bucket = read_bucket_list[bucket_idx];
            target_bucket.reserve(target_bucket.size() + source_bucket.size());
            std::ranges::move(source_bucket, std::back_inserter(target_bucket));
            ++bucket_idx;
        });
    });

    return std::pair{bin_size, std::move(read_bucket_list)};
}

} // namespace jstmap
