// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <libcontrib/seqan/alphabet.hpp>
#include <libcontrib/seqan/horspool_pattern.hpp>

TEST(horspool_pattern_test, dna4_pattern)
{
    using jst::contrib::operator""_dna4;
                   //0         1         2         3         4
                   //012345678901234567890123456789012345678901234
    auto haystack = "ACGTGACTAGCACGTGACTAGCACGTGACTAGCACGTGACTAGC"_dna4;
    auto needle = "GCACG"_dna4;

    static_assert(seqan3::cerealisable<jst::contrib::dna5>);

    jst::contrib::horspool_pattern pattern{needle};

    auto op = libjst::search_operation_old(pattern);

    std::vector<size_t> actual_positions{};
    op(std::views::all(haystack), [&] (auto const & finder) {
        actual_positions.push_back(seqan::beginPosition(finder));
    });
    EXPECT_RANGE_EQ(actual_positions, (std::vector{9, 20, 31}));
}
