// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief The method to create the index.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <jstmap/global/jstmap_types.hpp>
#include <jstmap/index/options.hpp>

namespace jstmap
{

seqan3::interleaved_bloom_filter<> create_index(rcs_store_t const &, index_options const &);

} // namespace jstmap
