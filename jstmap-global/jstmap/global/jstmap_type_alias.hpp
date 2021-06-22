// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides globally defined type aliases for the jstmap tools.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <libjst/journaled_sequence_tree.hpp>
#include <libjst/journal_sequence_tree_partitioned.hpp>

namespace jstmap
{

//!\brief The sequence type loaded from the disk.
using raw_sequence_t = std::vector<seqan3::dna5>;
//!\brief The globally available journal sequence tree type.
using jst_t = libjst::journaled_sequence_tree<raw_sequence_t>;
//!\brief The globally available partitioned journal sequence tree type.
using partitioned_jst_t = libjst::journal_sequence_tree_partitioned<jst_t>;

//!\brief An enum to select the verbosity level.
enum class verbosity_level : uint8_t
{
    quite, //!< No logging output is emitted.
    standard, //!< Logs regular information with no extra information on the output.
    verbose //!< Extra verbose logging output for debugging purposes.
};

//!\brief An enum to select the logging level.
enum class logging_level : uint8_t
{
    info, //!< An informative message during the execution.
    warning, //!< A warning message for non-severe issues during the execution.
    error //!< An error message for severe issues during the execution.
};

}  // namespace jstmap
