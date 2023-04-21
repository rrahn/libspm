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
#include <libjst/matcher/myers_matcher_restorable.hpp>

using jst::contrib::operator""_dna4;

struct myers_matcher_restorable_test : public ::testing::Test {
    using sequence_t = std::vector<jst::contrib::dna4>;
                         //0         1         2         3         4
                         //012345678901234567890123456789012345678901234
    sequence_t haystack = "ACGTGACTAGCACGTGACTAGCACGTGACTAGCACGTGACTAGC"_dna4;
    sequence_t needle = "GCACG"_dna4;
    std::size_t errors = 1;

    std::vector<std::size_t> expected_positions{13,14,15,24,25,26,35,36,37};

    auto get_matcher() const noexcept {
        return libjst::restorable_myers_matcher{needle, errors};
    }
};

TEST_F(myers_matcher_restorable_test, concept_tests) {
    using matcher_t = decltype(get_matcher());
    EXPECT_TRUE(libjst::window_matcher<matcher_t>);
}

TEST_F(myers_matcher_restorable_test, window_size) {
    auto matcher = get_matcher();
    EXPECT_EQ(libjst::window_size(matcher), std::ranges::size(needle) + errors);
}

TEST_F(myers_matcher_restorable_test, dna4_pattern)
{
    auto matcher = get_matcher();

    std::vector<size_t> actual_positions{};
    matcher(haystack, [&] (auto const & finder) {
        actual_positions.push_back(seqan::endPosition(finder));
    });
    EXPECT_RANGE_EQ(actual_positions, expected_positions);
}

TEST_F(myers_matcher_restorable_test, dna4_pattern_captured)
{
    size_t chunk_size{13};
    std::ptrdiff_t chunk_count = (haystack.size() + chunk_size - 1)/chunk_size;

    auto matcher = get_matcher();
    auto state = matcher.capture();
    std::vector<size_t> actual_positions{};
    for (std::ptrdiff_t chunk_idx = 0; chunk_idx < chunk_count; ++chunk_idx) {
        std::ptrdiff_t offset = chunk_idx * chunk_size;
        sequence_t chunk{haystack.begin() + offset,
                         haystack.begin() + std::min<std::ptrdiff_t>(offset + chunk_size, haystack.size())};
        matcher.restore(state);
        matcher(chunk, [&] (auto const & finder) {
            actual_positions.push_back(seqan::endPosition(finder) + offset);
        });
        state = matcher.capture();
    }
    EXPECT_RANGE_EQ(actual_positions, expected_positions);
}
