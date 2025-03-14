// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>

#include <libspm/seqan/alphabet.hpp>

#include <libspm/matcher/concept.hpp>
#include <libspm/matcher/myers_matcher.hpp>

using spm::operator""_dna4;

struct myers_matcher_test : public ::testing::Test {
    using sequence_t = std::vector<spm::dna4>;
                         //0         1         2         3         4
                         //012345678901234567890123456789012345678901234
    sequence_t haystack = "ACGTGACTAGCACGTGACTAGCACGTGACTAGCACGTGACTAGC"_dna4;
    sequence_t needle = "GCACG"_dna4;
    std::size_t errors = 1;

    std::vector<std::size_t> expected_positions{13,14,15,24,25,26,35,36,37};

    auto get_matcher() const noexcept {
        return spm::myers_matcher{needle, errors};
    }
};

TEST_F(myers_matcher_test, concept_tests) {
    using matcher_t = decltype(get_matcher());
    EXPECT_TRUE(spm::window_matcher<matcher_t>);
}

TEST_F(myers_matcher_test, window_size) {
    auto matcher = get_matcher();
    EXPECT_EQ(spm::window_size(matcher), std::ranges::size(needle) + errors);
}

TEST_F(myers_matcher_test, dna4_pattern)
{
    auto matcher = get_matcher();

    std::vector<size_t> actual_positions{};
    matcher(haystack, [&] (auto const & finder) {
        actual_positions.push_back(seqan2::endPosition(finder));
    });
    EXPECT_TRUE(std::ranges::equal(actual_positions, expected_positions));
}
