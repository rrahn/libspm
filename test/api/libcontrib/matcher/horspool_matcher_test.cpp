// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>

#include <libcontrib/seqan/alphabet.hpp>

#include <libcontrib/matcher/concept.hpp>
#include <libcontrib/matcher/horspool_matcher.hpp>

using jst::contrib::operator""_dna4;

struct horspool_matcher_test : public ::testing::Test {
    using sequence_t = std::vector<jst::contrib::dna4>;
                         //0         1         2         3         4
                         //012345678901234567890123456789012345678901234
    sequence_t haystack = "ACGTGACTAGCACGTGACTAGCACGTGACTAGCACGTGACTAGC"_dna4;
    sequence_t needle = "GCACG"_dna4;

    std::vector<std::size_t> expected_positions{9, 20, 31};

    auto get_matcher() const noexcept {
        return jst::contrib::horspool_matcher{needle};
    }
};

TEST_F(horspool_matcher_test, concept_tests) {
    using matcher_t = decltype(get_matcher());
    EXPECT_TRUE(jst::contrib::window_matcher<matcher_t>);
}

TEST_F(horspool_matcher_test, window_size) {
    auto matcher = get_matcher();
    EXPECT_EQ(jst::contrib::window_size(matcher), std::ranges::size(needle));
}

TEST_F(horspool_matcher_test, dna4_pattern)
{
    auto matcher = get_matcher();

    std::vector<size_t> actual_positions{};
    matcher(haystack, [&] (auto const & finder) {
        actual_positions.push_back(seqan2::beginPosition(finder));
    });
    EXPECT_TRUE(std::ranges::equal(actual_positions, expected_positions));
}

