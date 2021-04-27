// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>

#include <seqan3/test/expect_range_eq.hpp>

#include <jstmap/search/load_jst.hpp>

#include <libjst/journal_sequence_tree_partitioned.hpp>

TEST(journal_sequence_tree_partitioned_test, jst_partitioned)
{
    std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    jstmap::jst_t jst = jstmap::load_jst(jst_file);
    size_t bin_count = 3;
    libjst::journal_sequence_tree_partitioned jst_part(&jst, bin_count);

    EXPECT_EQ(jst_part.size(), bin_count);


    libjst::detail::journal_sequence_tree_context_enumerator bin_0_1 = jst_part[0];
    libjst::detail::journal_sequence_tree_context_enumerator bin_0_2 = jst_part[0];
    auto it1 = bin_0_1.begin();
    auto it2 = bin_0_2.begin();
    // EXPECT_EQ(bin_0_1, bin_0_2);
    for (size_t i = 0; i < 1; ++i)
    {
        ++it1;
        ++it2;
        // EXPECT_EQ(bin_0_1, bin_0_2);
        EXPECT_RANGE_EQ(*it1, *it2);
    }


    // size_t bin_size = 0;
    // for (auto it = jst_part[0].begin(); it != jst_part[0].end(); ++it)
    //     ++bin_size;
    //
    // EXPECT_EQ(bin_size, ceil(reference_size / bin_count))
}

TEST(journal_sequence_tree_partitioned_test, empty_jst_partitioned)
{

}

TEST(journal_sequence_tree_partitioned_test, jst_partitioned_negative_bin_count)
{

}
