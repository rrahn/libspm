// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <jstmap/search/load_jst.hpp>
#include <jstmap/search/load_queries.hpp>
#include <jstmap/search/search_queries.hpp>

TEST(jstmap_index, search_jst)
{
    std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    jstmap::jst_t jst = jstmap::load_jst(jst_file);

    std::filesystem::path queries_file{DATADIR"sim_reads_ref1x10.fa"};
    std::vector reads = jstmap::load_queries(queries_file);

    std::vector results = jstmap::search_queries(std::move(jst), std::move(reads));

    EXPECT_EQ(results[0], (libjst::context_position{0, 0}));
    EXPECT_EQ(results[1], (libjst::context_position{0, 16}));
    EXPECT_EQ(results[2], (libjst::context_position{0, 36}));
    EXPECT_EQ(results[3], (libjst::context_position{0, 1}));
    EXPECT_EQ(results[4], (libjst::context_position{0, 21}));
    EXPECT_EQ(results[5], (libjst::context_position{0, 41}));
    EXPECT_EQ(results[6], (libjst::context_position{0, 61}));
    EXPECT_EQ(results[7], (libjst::context_position{0, 70}));
    EXPECT_EQ(results[8], (libjst::context_position{0, 41}));
    EXPECT_EQ(results[9], (libjst::context_position{0, 50}));
}
