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
#include <seqan3/utility/range/to.hpp>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/all_matches.hpp>
#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/bam_writer.hpp>
#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/global/search_matches.hpp>
#include <jstmap/search/filter_queries.hpp>
#include <jstmap/search/match_aligner.hpp>
#include <jstmap/search/load_queries.hpp>
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

    search_parser.add_flag(options.is_quite,
                           'q',
                           "quite",
                           "Disables all logging.",
                           seqan3::option_spec::standard);
    search_parser.add_flag(options.is_verbose,
                           'v',
                           "verbose",
                           "Enables expansive debug logging.",
                           seqan3::option_spec::standard);

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
        if (options.is_quite) {
            get_application_logger().set_verbosity(verbosity_level::quite);
        } else if (options.is_verbose) {
            get_application_logger().set_verbosity(verbosity_level::verbose);
        }

        log_debug("References file:", options.jst_input_file_path.string());
        log_debug("Query file:", options.query_input_file_path.string());
        log_debug("Output file:", options.map_output_file_path.string());
        log_debug("Index file:", options.index_input_file_path.string());
        log_debug("Error rate:", options.error_rate);
        log_debug("Thread count:", options.thread_count);
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        log_err(ex.what());
        return -1;
    }

    // Load the queries and the jst.
    auto global_start = std::chrono::high_resolution_clock::now();
    try
    {
        log_info("Start mapping");
        log_debug("Load reads");
        auto start = std::chrono::high_resolution_clock::now();
        std::vector query_records = load_queries(options.query_input_file_path);
        size_t query_idx{};
        auto queries = query_records
                     | std::views::transform([&] (sequence_record_t & record) {
                        return search_query{query_idx++, std::move(record)};
                     })
                     | seqan3::ranges::to<std::vector>();
        auto end = std::chrono::high_resolution_clock::now();
        log_debug("Read count", queries.size());
        log_debug("Loading time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        log_debug("Load reference database");
        start = std::chrono::high_resolution_clock::now();
        rcs_store_t rcs_store = load_jst(options.jst_input_file_path);
        end = std::chrono::high_resolution_clock::now();
        log_debug("Loading time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        // Now we want to handle the bins as well.

        start = std::chrono::high_resolution_clock::now();
        size_t bin_size{std::numeric_limits<size_t>::max()};
        std::vector<bucket_type> bucket_list{};
        // What if the ibf is not present?
        if (options.index_input_file_path.empty())
        {
            log_debug("No prefilter enabled");
            bucket_list.resize(1);
            bucket_list[0] = queries;
        }
        else
        {
            log_debug("Applying IBF prefilter");
            std::tie(bin_size, bucket_list) = filter_queries(queries, options);
            log_debug("Bin size:", bin_size);
            log_debug("Bucket count:", bucket_list.size());
            // auto bucket_sizes = bucket_list | std::views::transform([](auto const & bucket) {
            //     return std::ranges::size(bucket);
            // });
            // auto non_empty_buckets = bucket_sizes | std::views::filter([](size_t const s) { return s > 0;});
            // log_debug("Non-empty bucket count:", std::ranges::distance(non_empty_buckets.begin(), non_empty_buckets.end()));
            // size_t candidates = std::accumulate(bucket_sizes.begin(), bucket_sizes.end(), 0);
            // log_debug("Candidate count:", candidates);
            // log_debug("Bucket sizes:", bucket_sizes);
        }
        end = std::chrono::high_resolution_clock::now();
        log_debug("Filter time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

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

        start = std::chrono::high_resolution_clock::now();
        auto query_matches = queries
                           | std::views::transform([] (search_query & query) {
                                return all_matches{std::move(query)};
                           })
                           | seqan3::ranges::to<std::vector>();


        for (size_t bucket_idx = 0; bucket_idx < bucket_list.size(); ++bucket_idx)
        { // parallel region
            auto const & bucket = bucket_list[bucket_idx];

            if (bucket.empty())
                continue;

            // Step 1: distribute search:
            log_debug("Local search in bucket:", bucket_idx);

            // Step 2: transform to haystack
            auto largest_query_it = std::ranges::max_element(bucket, std::ranges::less{}, [] (search_query const & query) {
                return std::ranges::size(query.value().sequence());
            });
            size_t const max_window_size = std::ranges::size((*largest_query_it).value().sequence()) - 1;
            auto chunked_tree = libjst::make_volatile(rcs_store) | libjst::chunk(bin_size, max_window_size);
            auto haystack = chunked_tree[bucket_idx];

            // Step 3: select matching strategy

            auto op = matching_operation{};
                // requires: options.error_count, bucket,

            // Step 4: apply matching
            op(std::move(haystack), bucket, [&] (search_query const & query, match_position position) {
                log_debug("Record match for query ", query.key(), " at ", position);
                query_matches[query.key()].record_match(std::move(position));
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
        log_debug("Matching time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");

        // Step 5: postprocess matches
        start = std::chrono::high_resolution_clock::now();

        std::vector<search_matches> aligned_matches_list{};
        aligned_matches_list.reserve(query_matches.size());

        std::ranges::for_each(query_matches, [&] (all_matches & query_match) {
            search_matches aligned_matches{std::move(query_match).query()};
            match_aligner aligner{rcs_store, aligned_matches.query().value().sequence()};

            for (match_position pos : query_match.matches()) {
                aligned_matches.record_match(aligner(std::move(pos)));
            }
            aligned_matches_list.push_back(std::move(aligned_matches));
        });

        end = std::chrono::high_resolution_clock::now();
        log_debug("Aligning time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
        // Filter globally for the best hits.
        // Ater reducing them to the same file sort, uniquify and then erase redundant matches.

        // Step 6: finalise
        start = std::chrono::high_resolution_clock::now();
        bam_writer writer{rcs_store, options.map_output_file_path};
        std::ranges::for_each(aligned_matches_list, [&] (search_matches const & matches) {
            writer.write_matches(matches);
        });
        end = std::chrono::high_resolution_clock::now();
        log_debug("Writing time:", std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), "s");
    }
    catch (std::exception const & ex)
    {
        log_err(ex.what());
        return -1;
    }
    auto global_end = std::chrono::high_resolution_clock::now();
    log_info("Finished mapping [", std::chrono::duration_cast<std::chrono::seconds>(global_end - global_start).count() ,"s]");
    return 0;
}

} // namespace jstmap
