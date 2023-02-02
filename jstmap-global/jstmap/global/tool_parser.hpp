// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides tool parser.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/argument_parser/argument_parser.hpp>

#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/options_base.hpp>

namespace jstmap
{

    inline void add_base_options(seqan3::argument_parser & parser, options_base & base_options) {
        parser.add_flag(base_options.is_quite,
                             'q',
                             "quite",
                             "Disables all logging.",
                             seqan3::option_spec::standard);
        parser.add_flag(base_options.is_verbose,
                             'v',
                             "verbose",
                             "Enables expansive debug logging.",
                             seqan3::option_spec::standard);
    }

    inline void initialise_logging_level(options_base const & base_options) noexcept {
        if (base_options.is_quite) {
            get_application_logger().set_verbosity(verbosity_level::quite);
        } else if (base_options.is_verbose) {
            get_application_logger().set_verbosity(verbosity_level::verbose);
        }
    }
}  // namespace jstmap
