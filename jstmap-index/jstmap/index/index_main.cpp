// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the main entry point of the just_map indexer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <seqan3/argument_parser/exceptions.hpp>

#include <jstmap/index/index_main.hpp>
#include <jstmap/index/options.hpp>

namespace jstmap
{

int index_main(seqan3::argument_parser & index_parser)
{
    index_options options{};

    index_parser.add_positional_option(options.input_file, "The input file.");
    index_parser.add_positional_option(options.output_file, "The output file.");

    try
    {
        index_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        return -1;
    }

    return 0;
}

} // namespace jstmap
