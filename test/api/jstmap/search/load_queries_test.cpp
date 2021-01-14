// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <jstmap/search/load_queries.hpp>

TEST(jstmap_index, load_queries)
{
    std::filesystem::path queries_file{DATADIR"sim_reads_ref1x10.fa"};
    std::vector reads = jstmap::load_queries(queries_file);

    EXPECT_EQ(reads.size(), 10u);
}
