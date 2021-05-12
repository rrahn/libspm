// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides globally defined types for the simulate subprogram.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <utility> //pair
#include <vector>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <libjst/journaled_sequence_tree.hpp>
#include <libjst/journal_sequence_tree_partitioned.hpp>

namespace jstmap
{

using sequence_t = std::vector<seqan3::dna5>;
using aligned_sequence_t = std::vector<seqan3::gapped<std::ranges::range_value_t<sequence_t>>>;
using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;
//!\brief The globally available journal sequence tree type.
using jst_t = libjst::journaled_sequence_tree<sequence_t>;
//!\brief The globally available partitioned journal sequence tree type.
using partitioned_jst_t = libjst::journal_sequence_tree_partitioned<jst_t>;

}  // namespace jstmap
