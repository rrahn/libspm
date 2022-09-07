// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/variant/alternate_sequence_kind.hpp>

struct alternate_sequence_kind_test : public ::testing::Test {
    static constexpr libjst::alternate_sequence_kind ins{libjst::alternate_sequence_kind::insertion};
    static constexpr libjst::alternate_sequence_kind rep{libjst::alternate_sequence_kind::replacement};
    static constexpr libjst::alternate_sequence_kind del{libjst::alternate_sequence_kind::deletion};
};

TEST_F(alternate_sequence_kind_test, insertion) {
    EXPECT_LT(ins, rep);
    EXPECT_LT(ins, del);
    EXPECT_EQ(ins, ins);
}

TEST_F(alternate_sequence_kind_test, replacement) {
    EXPECT_GT(rep, ins);
    EXPECT_LT(rep, del);
    EXPECT_EQ(rep, rep);
}

TEST_F(alternate_sequence_kind_test, deletion) {
    EXPECT_GT(del, ins);
    EXPECT_GT(del, rep);
    EXPECT_EQ(del, del);
}
