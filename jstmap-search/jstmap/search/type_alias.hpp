// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <ranges>

#include <seqan/sequence.h>

#include <libjst/sequence_tree/chunked_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>

#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/global/search_query.hpp>

namespace jstmap
{

using bin_sequence_t = std::views::all_t<reference_t const &>;
using bin_t = seqan::StringSet<bin_sequence_t>;

// Query types
using bucket_type = std::vector<search_query>;

// Haystack types
using jst_type = libjst::volatile_tree<rcs_store_t>;
using chunked_jst_type = decltype(std::declval<jst_type>() | libjst::chunk(1u));
using haystack_type = std::ranges::range_reference_t<chunked_jst_type>;

} // namespace jstmap
