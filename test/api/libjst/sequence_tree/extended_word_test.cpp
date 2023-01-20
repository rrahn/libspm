// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/sequence_tree/path_descriptor.hpp>

struct extended_word_test : public ::testing::Test {

};

TEST(extended_word_test, size) {
    libjst::extended_word<> word0{};
    EXPECT_EQ(word0.size(), 0ull);

    libjst::extended_word<> word7{0b001101010};
    EXPECT_EQ(word7.size(), 7ull);

    libjst::extended_word<> word70{word7};
    word70 <<= 63;
    EXPECT_EQ(word70.size(), 70ull);

    libjst::extended_word<> word133{word70};
    word133 <<= 63;
    EXPECT_EQ(word133.size(), 133ull);

    libjst::extended_word<> word196{word133};
    word196 <<= 63;
    EXPECT_EQ(word196.size(), 196ull);
}

TEST(extended_word_test, max_size) {
    libjst::extended_word<> word0{};
    EXPECT_EQ(word0.max_size(), 256ull);
}

TEST(extended_word_test, binary_or) {
    libjst::extended_word<> word{};
    word |= 1;
    EXPECT_EQ(word.size(), 1ull);
    EXPECT_EQ(word[0], true);
    word <<= 4;
    EXPECT_EQ(word.size(), 5ull);
    EXPECT_EQ(word[0], false);
    word |= 1;
    EXPECT_EQ(word.size(), 5ull);
    EXPECT_EQ(word[0], true);
    EXPECT_EQ(word[4], true);
}
