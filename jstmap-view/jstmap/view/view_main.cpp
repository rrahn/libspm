// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map viewer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

#include <jstmap/view/load_jst.hpp>
#include <jstmap/view/options.hpp>
#include <jstmap/view/view_format_fasta.hpp>

namespace jstmap
{

int view_main(seqan3::argument_parser & view_parser)
{
    view_options options{};

    view_parser.add_positional_option(options.jst_file,
                                       "The jst file.",
                                       seqan3::input_file_validator{{"jst"}});
    view_parser.add_option(options.haplotype_index,
                            '\0',
                            "haplotype-index",
                            "The index of the haplotype to print to the command line in fasta format.",
                            seqan3::option_spec::standard,
                            seqan3::arithmetic_range_validator<size_t>{0, std::numeric_limits<size_t>::max()});

    try
    {
        view_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ex)
    {
        std::cerr << "[ERROR] While parsing command line options: " << ex.what() << "\n";
        return -1;
    }

    // ----------------------------------------------------------------------------
    // Run the viewer.
    // ----------------------------------------------------------------------------

    int error_code = 0;
    try
    {
        auto jst = jstmap::load_jst(options.jst_file);
        view_as_format(jst, options.haplotype_index);
    }
    catch (std::exception const & ex)
    {
        std::cerr << "[ERROR] While viewing the jst content: " << ex.what() << "\n";
        error_code = -1;
    }

    return error_code;
}

} // namespace jstmap
