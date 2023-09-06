// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/matcher/concept.hpp>
#include <libjst/matcher/pigeonhole_matcher.hpp>

using jst::contrib::operator""_dna4;

struct pigeonhole_matcher_test : public ::testing::Test {
    using sequence_t = std::vector<jst::contrib::dna4>;
    using needle_position_t = seqan::PigeonholeSeedOnlyPosition;
                         //0         1         2         3         4
                         //012345678901234567890123456789012345678901234
    sequence_t haystack = "ACGTGACTAGCACGTGACTAGCACGTGACTAGCACGTGACTAGC"_dna4;
    sequence_t needle = "GCACG"_dna4;
    sequence_t needle2 = "TGACTAGCAC"_dna4;
    std::vector<sequence_t> multi_needle{needle, needle2};
    double errors = 0.0;

    std::vector<std::size_t> expected_positions{9, 20, 31};
    std::vector<std::size_t> expected_multi_positions{3, 8, 9, 14, 19, 20, 25, 30, 31, 36};
    std::vector<needle_position_t> expected_needle_positions{{1, 0, 5}, {1, 5, 5}, {0, 0, 5}, {1, 0, 5}, {1, 5, 5},
                                                             {0, 0, 5}, {1, 0, 5}, {1, 5, 5}, {0, 0, 5}, {1, 0, 5}};

    auto get_matcher() const noexcept {
        return libjst::pigeonhole_matcher{needle, errors};
    }

    auto get_multi_matcher() const noexcept {
        return libjst::pigeonhole_matcher{multi_needle, errors};
    }
};

TEST_F(pigeonhole_matcher_test, concept_tests) {
    using matcher_t = decltype(get_matcher());
    EXPECT_TRUE(libjst::window_matcher<matcher_t>);
}

TEST_F(pigeonhole_matcher_test, window_size) {
    auto matcher = get_matcher();
    EXPECT_EQ(libjst::window_size(matcher), std::ranges::size(needle));
}

TEST_F(pigeonhole_matcher_test, dna4_pattern)
{
    auto matcher = get_matcher();

    std::vector<size_t> actual_positions{};
    matcher(haystack, [&] (auto const & finder) {
        actual_positions.push_back(seqan::beginPosition(finder));
    });
    EXPECT_RANGE_EQ(actual_positions, expected_positions);
}

TEST_F(pigeonhole_matcher_test, dna4_multi_pattern)
{
    auto matcher = get_multi_matcher();

    std::vector<size_t> actual_positions{};
    matcher(haystack, [&] (auto const & finder) {
        actual_positions.push_back(seqan::beginPosition(finder));
    });
    EXPECT_RANGE_EQ(actual_positions, expected_multi_positions);
}

TEST_F(pigeonhole_matcher_test, dna4_multi_pattern_position)
{
    auto matcher = get_multi_matcher();

    std::vector<needle_position_t> actual_needle_positions{};
    matcher(haystack, [&] (auto const &) {
        actual_needle_positions.push_back(matcher.position());
    });
    EXPECT_EQ(actual_needle_positions, expected_needle_positions);
}
