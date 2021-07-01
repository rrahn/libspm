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
#include <omp.h>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <jstmap/global/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/filter_queries.hpp>
#include <jstmap/search/search_main.hpp>
#include <jstmap/search/search_queries.hpp>
#include <jstmap/search/type_alias.hpp>
#include <jstmap/search/write_results.hpp>
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
        auto queries = load_queries(options.query_input_file_path);
        auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Load queries time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        std::cout << "load the jst\n";
        start = std::chrono::high_resolution_clock::now();
        auto jst = load_jst(options.jst_input_file_path);
        end = std::chrono::high_resolution_clock::now();

        std::cout << "Load JST time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";

        // Now we want to handle the bins as well.
        // Step 1: prepare bins: seqan::StringSet<std::views::all_t<raw_sequence_t &>>
        // Step 2:
        start = std::chrono::high_resolution_clock::now();
        size_t bin_size{std::numeric_limits<size_t>::max()};
        std::vector<bin_t> bins{};
        // What if the ibf is not present?
        if (options.index_input_file_path.empty())
        {
            bins.resize(1);
            std::ranges::for_each(queries, [&] (raw_sequence_t const & query)
            {
                seqan::appendValue(bins.front(), query | std::views::all);
            });
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

        std::cout << "bin_size = " << bin_size << "\n";
        start = std::chrono::high_resolution_clock::now();
        partitioned_jst_t pjst{std::addressof(jst), bin_size}; // default initialisation
        std::cout << "bin count = " << pjst.bin_count() << "\n";
        std::cout << "Run search\n";

        // * filter step with ibf -> {bin_id, {ref_view(query_l)[, ref_view(query_r)], global_query_id}[]}
        // list of {bin_id:queries}
        // partioned_jst[bin_id] -> traverser_model:
            // range_agent{traverser_model, } we can construct this from the model directly.

        // We need to write more information including the reference sequences and the length of the reference sequences.
        std::vector<std::vector<search_match>> bin_matches{};
        bin_matches.resize(pjst.bin_count());

        #pragma omp parallel for num_threads(options.thread_count) shared(bin_matches, pjst, bins, options) schedule(dynamic)
        for (size_t bin_idx = 0; bin_idx < pjst.bin_count(); ++bin_idx)
        { // parallel region
            if (seqan::empty(bins[bin_idx]))
                continue;

            auto jst_bin = pjst.bin_at(bin_idx);
            // * search queries in bin_id -> matches[]
            // * push results into global queue
            // seqan::StringSet<raw_sequence_t> _queries{};
            // seqan::appendValue(_queries, queries.front());
            // seqan::reserve(_queries, queries.size());

            // for (unsigned i = 0; i < queries.size(); ++i)
            //     seqan::appendValue(_queries, queries[i]);
            // std::cout << "Searching bin = " << bin_idx << "\n";
            bin_matches[bin_idx] = search_queries_(jst_bin, bins[bin_idx], options.error_rate);

            std::ranges::sort(bin_matches[bin_idx], std::less<void>{}); // sort by query_id and by error_count.
            auto redundant_tail = std::ranges::unique(bin_matches[bin_idx],
                                                      std::ranges::equal_to{},
                                                      [] (auto const & match) { return match.query_id; });
            bin_matches[bin_idx].erase(redundant_tail.begin(), redundant_tail.end());

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

        std::cout << "Write results\n";
        start = std::chrono::high_resolution_clock::now();
        sam_file_t sam_file{options.map_output_file_path};

        for (size_t bin_idx = 0; bin_idx < pjst.bin_count(); ++bin_idx)
            if (!bin_matches[bin_idx].empty())
                write_results(sam_file, bin_matches[bin_idx], bins[bin_idx]); // needs synchronised output_buffer

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Write results time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  << " sec\n";
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
