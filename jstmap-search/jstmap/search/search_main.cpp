// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map searcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <functional>
#include <numeric>
#include <omp.h>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/filter_queries.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/match.hpp>
#include <jstmap/search/matching_operation.hpp>
#include <jstmap/search/search_main.hpp>
// #include <jstmap/search/search_queries.hpp>
#include <jstmap/search/type_alias.hpp>
// #include <jstmap/search/write_results.hpp>
#include <jstmap/search/options.hpp>

namespace jstmap
{

int search_main(seqan3::argument_parser & search_parser)
{
    search_options options{};

    search_parser.add_positional_option(options.jst_input_file_path,
                                       "The path to the journaled sequence tree.",
                                       seqan3::input_file_validator{{"jst"}});
    search_parser.add_positional_option(options.query_input_file_path,
                                       "The path to the read file.",
                                       seqan3::input_file_validator{{"fa", "fasta"}});
    search_parser.add_positional_option(options.map_output_file_path,
                                       "The alignment map output file.",
                                       seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                     {"sam", "bam"}});
    search_parser.add_option(options.index_input_file_path,
                             'i',
                             "index",
                             "The prebuilt index to speedup the search.",
                             seqan3::option_spec::standard,
                             seqan3::input_file_validator{{"ibf"}});
    search_parser.add_option(options.error_rate,
                             'e',
                             "error-rate",
                             "The error rate allowed for mapping the reads.",
                             seqan3::option_spec::standard,
                             seqan3::arithmetic_range_validator{0.0, 1.0});
    search_parser.add_option(options.thread_count,
                             't',
                             "thread-count",
                             "The number of threads to use for the search.",
                             seqan3::option_spec::standard,
                             seqan3::arithmetic_range_validator{1u, std::thread::hardware_concurrency()});

    try
    {
        search_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    // Load the queries and the jst.
    try
    {
        std::cout << "load the queries\n";
        auto start = std::chrono::high_resolution_clock::now();
        std::vector const queries = load_queries(options.query_input_file_path);
        auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Load queries time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        std::cout << "load the jst\n";
        start = std::chrono::high_resolution_clock::now();

        rcs_store_t rcs_store = load_jst(options.jst_input_file_path);
        // seqan3::debug_stream << "Reference: " << rcs_store.source() << "\n";
        seqan3::debug_stream << "Size: " << rcs_store.size() << "\n";

        end = std::chrono::high_resolution_clock::now();

        std::cout << "Load JST time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // Now we want to handle the bins as well.
        // Step 1: prepare bins: seqan::StringSet<std::views::all_t<raw_sequence_t &>>
        // Step 2:
        start = std::chrono::high_resolution_clock::now();
        size_t bin_size{std::numeric_limits<size_t>::max()};
        std::vector<bucket_type> bucket_list{};
        // What if the ibf is not present?
        if (options.index_input_file_path.empty())
        {
            bucket_list.resize(1);
            bucket_list[0].reserve(queries.size());
            query_index_type idx{};
            std::ranges::for_each(queries, [&] (auto const & query) {
                bucket_list[0].emplace_back(idx++, query | std::views::all);
                // seqan::appendValue(bucket_list[0], (decltype(query)&&)query | std::views::all);
            });
        }
        else
        {
            std::cout << "Filter the queries\n";
            std::tie(bin_size, bucket_list) = filter_queries(queries, options);
        }
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Filter time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // #TODO: add later
        // std::cout << "bin_size = " << bin_size << "\n";
        // start = std::chrono::high_resolution_clock::now();
        // partitioned_jst_t pjst{std::addressof(jst), bin_size}; // default initialisation
        // std::cout << "bin count = " << pjst.bin_count() << "\n";
        // std::cout << "Run search\n";

        // * filter step with ibf -> {bin_id, {ref_view(query_l)[, ref_view(query_r)], global_query_id}[]}
        // list of {bin_id:queries}
        // partioned_jst[bin_id] -> traverser_model:
            // range_agent{traverser_model, } we can construct this from the model directly.

        // We need to write more information including the reference sequences and the length of the reference sequences.
        std::vector<std::vector<match>> query_matches{};
        query_matches.resize(queries.size());

        auto chunked_tree = libjst::volatile_tree{rcs_store} | libjst::chunk(bin_size);

        for (size_t bucket_idx = 0; bucket_idx < bucket_list.size(); ++bucket_idx)
        { // parallel region
            auto const & bucket = bucket_list[bucket_idx];

            if (bucket.empty())
                continue;

            // Step 1: distribute search:

            // Step 2: transform to haystack
            auto haystack = chunked_tree[bucket_idx];

            // Step 3: select matching strategy

            auto op = matching_operation{};
                // requires: options.error_count, bucket,

            // Step 4: apply matching
            op(std::move(haystack), bucket, [&] (auto match) {
                query_matches[match.query_id()].push_back(std::move(match));
            });

        // }

        // * filter out duplicates? parallel?
        // * A) same read/read_pair is found in multiple locations (report only one hit per bin?)
        //  * depends on mapping mode: if best or all best then only the mapping locations with the lowest error count
        //  * if all then all alternative mapping locations sorted by their error count
        //      * needs: error_count and query_id for filtering
        // * B) If same hit identified in two bins because of bin-overlap.
        //  * depends on the HIBF setting: if with overlap then yes if not then search unidentified reads in overlap region with tight window.
        //  * how many reads are left and how many regions must be searched?
        //      * assume same base coordinate to filter

        // { // parallel region
            // * run simd alignment to obtain CIGAR string on all filtered matches and report
            //  -> multi-threaded conversion to record and synchronised buffer
                // -> one buffer per thread to fill -> full buffer is pushed into queue -> empty buffer is reserved when available
                // -> single buffer writes out record stream into bgzf_ostream (possibly unsynchronised?)
            // * report in BAM file
            // We do not want to have a critical region here.
            // But somehow we need to provide some mechanism to say there is a new batch available.
            // We need the concurrent queue here!

        }

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Search time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";
        size_t total_occurrences{};
        std::ranges::for_each(query_matches, [&] (auto const & matches) { total_occurrences += matches.size(); });
        std::cout << "Total hits: " << total_occurrences << "\n";
        std::cout << "Reduce and filter matches\n";
        start = std::chrono::high_resolution_clock::now();
        // Filter globally for the best hits.
        // Ater reducing them to the same file sort, uniquify and then erase redundant matches.

        // Step 5: postprocess matches

        // #TODO: readd again
        // std::vector<search_match2> total_matches{};
        // if (bucket_list.size() == 1)
        // {
        //     total_matches = std::move(bin_matches[0]);
        // }
        // else // Reduce all matches and uniquify all matches.
        // {
        //     auto match_counts = bin_matches | std::views::transform([] (auto & bin) { return bin.size(); });
        //     size_t total_match_count = std::accumulate(match_counts.begin(), match_counts.end(), 0);
        //     total_matches.resize(total_match_count);
        //     auto insert_it = total_matches.begin();
        //     std::ranges::for_each(bin_matches, [&] (auto && bin)
        //     {
        //         insert_it = std::ranges::move(bin, insert_it).out;
        //     });

        //     std::ranges::sort(total_matches, std::less<void>{}); // sort by query_id and by error_count.
        //     auto redundant_tail = std::ranges::unique(total_matches,
        //                                               std::ranges::equal_to{},
        //                                               [] (auto const & match) { return match.query_id; });
        //     total_matches.erase(redundant_tail.begin(), redundant_tail.end());
        // }

        end = std::chrono::high_resolution_clock::now();
        // std::cout << "Found " << total_matches.size() << " matches\n";
        std::cout << "Reduction time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // Step 6: finalise

        // std::cout << "Write results\n";
        // start = std::chrono::high_resolution_clock::now();
        // write_results(total_matches, queries, options); // needs synchronised output_buffer

        // end = std::chrono::high_resolution_clock::now();
        // std::cout << "Write results time: "
        //           << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        //           << " sec\n";

        // Step 7: write out
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
