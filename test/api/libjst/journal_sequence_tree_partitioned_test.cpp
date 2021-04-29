// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <libjst/journal_sequence_tree_partitioned.hpp>

#include "journal_sequence_tree_traversal_test_template.hpp"

struct partitioned_traversal_test : libjst::test::traversal_fixture_base
{};

TEST_P(partitioned_traversal_test, construct)
{
    auto jst = this->construct_jst();

    EXPECT_EQ(jst.size(), this->sequences.size());

    for (size_t i = 0; i < jst.size(); ++i)
        EXPECT_RANGE_EQ(jst.sequence_at(i), this->sequences[i]);
}

TEST_P(partitioned_traversal_test, enumerate_contexts)
{
    auto jst = this->construct_jst();
    libjst::journal_sequence_tree_partitioned p_jst{std::addressof(jst), GetParam().bin_count};

    auto context_enumerator = p_jst.context_enumerator(libjst::context_size{GetParam().context_size},
                                                       libjst::bin_index{0u});
    for (auto context_it = context_enumerator.begin(); context_it != context_enumerator.end(); ++context_it)
    {
        auto context = *context_it;
        std::string tmp = libjst::test::sequence_to_string(context);

        auto positions = context_it.positions();

        EXPECT_TRUE((this->context_positions_exist(tmp, positions))) << "context " << tmp;
    }

    // Verify that all unique contexts have been enumerated and that there is no unknown location.
    EXPECT_TRUE(this->all_contexts_enumerated());
    print_unvisited_contexts();

    EXPECT_TRUE(this->unknown_locations.empty());
    print_unknown_context_locations();
}

// ----------------------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(substitution_1, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //          0123456
    //               b
    // 0:       aaaa     [0, 0, 0, 0]
    // 1:        aaaa    [1, 1, 1, 1]
    // 2:         aaab   [-, 2, 2, -]
    // 3:          aaba  [-, 3, 3, -]
    // 4:         aaaa   [2, -, -, 2]
    // 5:          aaaa  [3, -, -, 3]
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{4u},
    .bin_count{1u}
}));
