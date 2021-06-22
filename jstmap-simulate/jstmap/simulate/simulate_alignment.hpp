// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <ranges>

#include <seqan3/alphabet/gap/gapped.hpp>

#include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap
{

//!\brief The aligned sequence type.
using aligned_sequence_t = std::vector<seqan3::gapped<std::ranges::range_value_t<raw_sequence_t>>>;
//!\brief The pairwise alignment type.
using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;

alignment_t simulate_alignment(raw_sequence_t & reference, double error_rate);

}
