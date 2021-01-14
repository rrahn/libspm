// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>

#include <jstmap/search/load_jst.hpp>

TEST(jstmap_index, load_jst)
{
    std::filesystem::path jst_file{DATADIR"sim_refx5.jst"};
    jstmap::jst_t jst = jstmap::load_jst(jst_file);

    EXPECT_EQ(jst.size(), 5u);
}

TEST(jstmap_index, load_jst_empty_path)
{
    EXPECT_THROW(jstmap::load_jst(std::filesystem::path{}), std::runtime_error);
}

TEST(jstmap_index, load_jst_unkown_path)
{
    std::filesystem::path unknown_jst_file{DATADIR"unknown.jst"};
    EXPECT_THROW(jstmap::load_jst(unknown_jst_file), std::runtime_error);
}
