// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides build function to create the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <jstmap/index/global_types.hpp>

namespace jstmap
{

//!\cond
std::pair<jst_t, partitioned_jst_t> build_journaled_sequence_tree(std::vector<raw_sequence_t> &&, 
                                                                            const uint32_t = 1);
//!\endcond

}  // namespace jstmap
