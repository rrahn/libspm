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

#include <jstmap/global/jstmap_jst_types.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/filter_queries.hpp>
#include <jstmap/search/search_main.hpp>
#include <jstmap/search/search_queries.hpp>
#include <jstmap/search/type_alias.hpp>
#include <jstmap/search/write_results.hpp>
#include <jstmap/search/options.hpp>

#include <libjst/journaled_sequence_tree/serialiser_concept.hpp>
#include <libjst/journaled_sequence_tree/serialiser_direct.hpp>
#include <libjst/journaled_sequence_tree/serialiser_delegate.hpp>

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
        // TODO replace for now!
        raw_sequence_t reference{};
        jst_model_t jst_model{reference, 0};
        fwd_jst_t jst{jst_model};
        // load the jst!
        {
            std::ifstream archive_stream{options.jst_input_file_path, std::ios_base::binary | std::ios_base::in};
            cereal::BinaryInputArchive in_archive{archive_stream};
            auto jst_archive = in_archive | libjst::direct_serialiser(reference)
                                        | libjst::delegate_serialiser(jst_model);
            libjst::load(jst, jst_archive);
        }
        end = std::chrono::high_resolution_clock::now();

        std::cout << "Load JST time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // Now we want to handle the bins as well.
        // Step 1: prepare bins: seqan::StringSet<std::views::all_t<raw_sequence_t &>>
        // Step 2:
        start = std::chrono::high_resolution_clock::now();
        size_t bin_size{std::numeric_limits<size_t>::max()};
        std::vector<std::vector<size_t>> bins{};
        // What if the ibf is not present?
        if (options.index_input_file_path.empty())
        {
            bins.resize(1);
            bins[0].resize(queries.size());
            std::iota(bins[0].begin(), bins[0].end(), 0ull);
        }
        else
        {
            std::cout << "Filter the queries\n";
            std::tie(bin_size, bins) = filter_queries(queries, options);
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
        std::vector<std::vector<search_match2>> bin_matches{};
        bin_matches.resize(1 /*pjst.bin_count()*/);

        #pragma omp parallel for num_threads(options.thread_count) shared(bin_matches, jst, bins, options) schedule(dynamic)
        for (size_t bin_idx = 0; bin_idx < 1; ++bin_idx)
        { // parallel region
            auto const & bin_query_ids = bins[bin_idx];

            if (bin_query_ids.empty())
                continue;

            // Step1: generate the local sequence set.
            bin_t local_bin{};
            seqan::reserve(local_bin, bin_query_ids.size());

            for (size_t local_bin_idx = 0; local_bin_idx < bin_query_ids.size(); ++local_bin_idx)
            {
                seqan::appendValue(local_bin, queries[bin_query_ids[local_bin_idx]] | std::views::all);
            }

            auto & jst_bin = jst; //pjst.bin_at(bin_idx);
            // * search queries in bin_id -> matches[]
            // * push results into global queue
            // seqan::StringSet<raw_sequence_t> _queries{};
            // seqan::appendValue(_queries, queries.front());
            // seqan::reserve(_queries, queries.size());

            // for (unsigned i = 0; i < queries.size(); ++i)
            //     seqan::appendValue(_queries, queries[i]);
            // std::cout << "Searching bin = " << bin_idx << "\n";
            bin_matches[bin_idx] = search_queries_horsppol(jst_bin, local_bin, options.error_rate);

            std::ranges::sort(bin_matches[bin_idx], std::less<void>{}); // sort by query_id and by error_count.
            auto redundant_tail = std::ranges::unique(bin_matches[bin_idx],
                                                      std::ranges::equal_to{},
                                                      [] (auto const & match) { return match.query_id; });
            bin_matches[bin_idx].erase(redundant_tail.begin(), redundant_tail.end());

            // Overwrite the query id here with the global query_id.
            std::ranges::for_each(bin_matches[bin_idx], [&] (search_match2 & match)
            {
                match.query_id = bins[bin_idx][match.query_id];
            });

            // seqan3::debug_stream << "Report " << matches.size() << " matches:\n";
            // for (auto const & match : matches)
            //     seqan3::debug_stream << "\t- match with: " << match.error_count << " errors and sequence: "
            //                          << match.sequence() << "\n";

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

        std::cout << "Reduce and filter matches\n";
        start = std::chrono::high_resolution_clock::now();
        // Filter globally for the best hits.
        // Ater reducing them to the same file sort, uniquify and then erase redundant matches.

        // #TODO: readd again
        std::vector<search_match2> total_matches{};
        if (bins.size() == 1)
        {
            total_matches = std::move(bin_matches[0]);
        }
        else // Reduce all matches and uniquify all matches.
        {
            auto match_counts = bin_matches | std::views::transform([] (auto & bin) { return bin.size(); });
            size_t total_match_count = std::accumulate(match_counts.begin(), match_counts.end(), 0);
            total_matches.resize(total_match_count);
            auto insert_it = total_matches.begin();
            std::ranges::for_each(bin_matches, [&] (auto && bin)
            {
                insert_it = std::ranges::move(bin, insert_it).out;
            });

            std::ranges::sort(total_matches, std::less<void>{}); // sort by query_id and by error_count.
            auto redundant_tail = std::ranges::unique(total_matches,
                                                      std::ranges::equal_to{},
                                                      [] (auto const & match) { return match.query_id; });
            total_matches.erase(redundant_tail.begin(), redundant_tail.end());
        }

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Found " << total_matches.size() << " matches\n";
        std::cout << "Reduction time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // std::cout << "Write results\n";
        // start = std::chrono::high_resolution_clock::now();
        // write_results(total_matches, queries, options); // needs synchronised output_buffer

        // end = std::chrono::high_resolution_clock::now();
        // std::cout << "Write results time: "
        //           << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        //           << " sec\n";
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
