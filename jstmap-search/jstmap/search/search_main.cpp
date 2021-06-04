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

#include <filesystem>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/search/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_main.hpp>
#include <jstmap/search/search_queries.hpp>
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
        auto queries = load_queries(options.query_input_file_path);

        std::cout << "load the jst\n";
        auto [jst, partitioned_jst] = load_jst(options.jst_input_file_path);

        std::vector results = search_queries(std::move(partitioned_jst), std::move(queries));

        write_results(std::move(results), options.map_output_file_path);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << "\n";
        return -1;
    }

    return 0;
}

} // namespace jstmap
