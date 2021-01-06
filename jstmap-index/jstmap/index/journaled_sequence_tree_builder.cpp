// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <jstmap/index/journaled_sequence_tree_builder.hpp>

namespace jstmap
{

jst_t build_journaled_sequence_tree(std::vector<raw_sequence_t> && sequences)
{
    return jst_t{std::move(sequences)};
}

}  // namespace jstmap
