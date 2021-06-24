// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <jstmap/create/journaled_sequence_tree_builder.hpp>

#include "test_data.hpp"

TEST(jstmap_index, build_partitioned_jst)
{
    std::vector<jstmap::raw_sequence_t> data{jstmap::test::sequences.size() + 1};
    data[0] = jstmap::raw_sequence_t{jstmap::test::reference};
    size_t index = 0;
    std::for_each(std::next(data.begin()), data.end(), [&] (jstmap::raw_sequence_t & sequence)
    {
        sequence = jstmap::test::sequences[index++];
    });

    auto [jst, partitioned_jst] = jstmap::build_journaled_sequence_tree(std::move(data), 2u);
    EXPECT_EQ(jst.size(), jstmap::test::sequences.size() + 1);
}
