// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/alphabet.hpp>
#include <libjst/rcms/delta_sequence_variant.hpp>

using jst::contrib::operator""_dna4;

struct delta_sequence_variant_test : public ::testing::Test {
    using source_t = std::vector<jst::contrib::dna4>;
    using sequence_type = libjst::delta_sequence_variant<source_t>;

};

TEST_F(delta_sequence_variant_test, concept) {
    EXPECT_TRUE(std::ranges::contiguous_range<sequence_type>);
    EXPECT_TRUE(std::ranges::view<sequence_type>);
    EXPECT_TRUE(std::ranges::sized_range<sequence_type>);
    EXPECT_TRUE(std::ranges::common_range<sequence_type>);
    EXPECT_TRUE(std::ranges::borrowed_range<sequence_type>);
}

TEST_F(delta_sequence_variant_test, snv) {
    jst::contrib::dna4 snv{'A'};
    sequence_type test{snv};
    EXPECT_RANGE_EQ(test, "A"_dna4);
}

TEST_F(delta_sequence_variant_test, deletion) {
    sequence_type test{};
    EXPECT_RANGE_EQ(test, source_t{});
}

TEST_F(delta_sequence_variant_test, insertion) {
    source_t src{"CGGACG"_dna4};
    sequence_type test{src};
    EXPECT_RANGE_EQ(test, src);
}
