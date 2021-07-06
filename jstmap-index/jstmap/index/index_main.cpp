// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map indexer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/load_jst.hpp>
#include <jstmap/index/create_index.hpp>
#include <jstmap/index/index_main.hpp>
#include <jstmap/index/options.hpp>
#include <jstmap/index/save_index.hpp>

namespace jstmap
{

int index_main(seqan3::argument_parser & index_parser)
{
    index_options options{};

    index_parser.add_positional_option(options.jst_input_file,
                                       "The input file.",
                                       seqan3::input_file_validator{{"jst"}});
    index_parser.add_positional_option(options.output_file,
                                       "The output file.",
                                       seqan3::output_file_validator{seqan3::output_file_open_options::create_new,
                                                                     {"ibf"}});
    index_parser.add_flag(options.is_quite,
                          '\0',
                          "quite",
                          "No logging output will be emitted.");
    index_parser.add_flag(options.is_verbose,
                          '\0',
                          "verbose",
                          "Verbose logging output will be emitted.");
    index_parser.add_option(options.bin_size,
                            'b',
                            "bin-size",
                            "The size of each bin for the index construction.",
                            seqan3::option_spec::standard);
    index_parser.add_option(options.kmer_size,
                            'k',
                            "kmer-size",
                            "The kmer-size used for the ibf creation.",
                            seqan3::option_spec::advanced,
                            seqan3::arithmetic_range_validator{0u, 31u});

    try
    {
        index_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        get_application_logger()(verbosity_level::standard, logging_level::error,
                                 "While parsing command line options: ", ex.what());
        return -1;
    }

    // ----------------------------------------------------------------------------
    // Initialise the global logger
    // ----------------------------------------------------------------------------

    jstmap::application_logger logger{true, options.is_quite
                                                ? verbosity_level::quite
                                                : options.is_verbose
                                                    ? verbosity_level::verbose
                                                    : verbosity_level::standard};
    set_application_logger(&logger);
    jstmap::application_logger & log = get_application_logger();

    // ----------------------------------------------------------------------------
    // Run the index creation.
    // ----------------------------------------------------------------------------

    log(verbosity_level::standard, logging_level::info, "Start index creation");

    int error_code = 0;
    try
    {
        log(verbosity_level::standard, logging_level::info, "Load jst: ", options.jst_input_file);
        auto jst = load_jst(options.jst_input_file);

        log(verbosity_level::standard, logging_level::info, "Creating the index with bin size ", options.bin_size,
                                                            " and kmer-size ", options.kmer_size);
        auto ibf = create_index(jst, options);

        log(verbosity_level::standard, logging_level::info, "Saving index: ", options.output_file);
        save_index(ibf, options);
    }
    catch (std::exception const & ex)
    {
        log(verbosity_level::standard, logging_level::error, "While creating the index: ", ex.what());
        error_code = -1;
    }
    log(verbosity_level::standard, logging_level::info, "Stop index creation");
    return error_code;
}

} // namespace jstmap
