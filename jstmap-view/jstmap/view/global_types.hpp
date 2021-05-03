// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides globally defined types for the view subprogram.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <libjst/journaled_sequence_tree.hpp>

namespace jstmap
{

//!\brief The sequence type loaded from the disk.
using raw_sequence_t = std::vector<seqan3::dna5>;
//!\brief The globally available journal sequence tree type.
using jst_t = libjst::journaled_sequence_tree<raw_sequence_t>;

}  // namespace jstmap
