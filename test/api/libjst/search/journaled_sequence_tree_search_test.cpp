// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/horspool_pattern.hpp>

#include <libjst/concept.hpp>
// #include <libjst/set.hpp>
#include <libjst/structure/jst_forward.hpp>

#include "../journal_sequence_tree_traversal_test_template.hpp"

struct jst_forward_test : public libjst::test::traversal_fixture_base
{};

// must be some type of receiver!
struct evaluate
{
    std::string _needle{};
    std::vector<size_t> _expected_positions{};
    std::vector<size_t> _actual_positions{};

    template <typename finder_t>
    void set_next(finder_t && finder) noexcept
    {
        _actual_positions.push_back(seqan::beginPosition(finder));
        EXPECT_RANGE_EQ(seqan::infix(finder), _needle); // in case of exact match!
    }

    void set_value() noexcept
    {
        EXPECT_EQ(_actual_positions.size(), _expected_positions.size());
        std::ranges::sort(_expected_positions);
        std::ranges::sort(_actual_positions);
        EXPECT_RANGE_EQ(_actual_positions, _expected_positions);
    }
};

TEST_P(jst_forward_test, construct)
{
    auto jst = this->construct_jst();

    EXPECT_EQ(jst.size(), this->sequences.size());

    for (size_t i = 0; i < jst.size(); ++i)
        EXPECT_RANGE_EQ(jst.sequence_at(i), this->sequences[i]);
}

TEST_P(jst_forward_test, search_horspool)
{
    auto jst = this->construct_jst();
    // jst.print_event_queue();
    libjst::journaled_sequence_tree_forward fwd_jst{std::move(jst)};

    jst::contrib::horspool_pattern pattern{GetParam().needle};
    auto sender = fwd_jst.search(pattern);
    auto op = sender.connect(evaluate{._needle = GetParam().needle,
                                      ._expected_positions = GetParam().expected_positions});
    op.start();
}

// ----------------------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(no_variants, jst_forward_test, testing::Values(
libjst::test::traversal_fixture
{
    //          0123456
    // 0:       aaaa     [0, 0, 0, 0]
    // 1:        aaaa    [1, 1, 1, 1]
    // 2:         aaaa   [2, 2, 2, 2]
    // 3:          aaaa  [3, 3, 3, 3]
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events{},
    .context_size{4}, // TODO: determine automatically
    .needle{"aaaa"s},
    .expected_positions{0, 1, 2, 3}
}));

INSTANTIATE_TEST_SUITE_P(single_snp_at_5, jst_forward_test, testing::Values(
libjst::test::traversal_fixture
{
    //             ____
    //          0123456
    //               b
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{4}, // TODO: determine automatically
    .needle{"aaba"s},
    .expected_positions{3}
}));

INSTANTIATE_TEST_SUITE_P(single_snp_1, jst_forward_test, testing::Values(
libjst::test::traversal_fixture
{
    //             ____
    //          0123456
    //               b
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{libjst::test::position_t{.offset = 5u}, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{4}, // TODO: determine automatically
    .needle{"aaaa"s},
    .expected_positions{0, 1, 2, 3}
}));

