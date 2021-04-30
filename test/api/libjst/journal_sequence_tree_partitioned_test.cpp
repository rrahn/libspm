// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <libjst/journal_sequence_tree_partitioned.hpp>

#include <cereal/archives/binary.hpp>

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

    EXPECT_GT(GetParam().bin_count, 0u);
    libjst::journal_sequence_tree_partitioned p_jst{std::addressof(jst), GetParam().bin_count};

    for (uint32_t index = 0; index < p_jst.bin_count(); ++index)
    {
        auto context_enumerator = p_jst.context_enumerator(libjst::context_size{GetParam().context_size},
                                                           libjst::bin_index{index});
        for (auto context_it = context_enumerator.begin(); context_it != context_enumerator.end(); ++context_it)
        {
            auto context = *context_it;
            std::string tmp = libjst::test::sequence_to_string(context);

            auto positions = context_it.positions();

            EXPECT_TRUE((this->context_positions_exist(tmp, positions))) << "context " << tmp;
        }
    }

    // Verify that all unique contexts have been enumerated and that there is no unknown location.
    EXPECT_TRUE(this->all_contexts_enumerated());
    print_unvisited_contexts();

    EXPECT_TRUE(this->unknown_locations.empty());
    print_unknown_context_locations();
}

TEST_P(partitioned_traversal_test, serialisation_test)
{
    auto jst = this->construct_jst();

    EXPECT_GT(GetParam().bin_count, 0u);

    libjst::journal_sequence_tree_partitioned p_jst_original{std::addressof(jst), GetParam().bin_count};
    // sereialise
    std::stringstream archive_stream{};
    {
        cereal::BinaryOutputArchive archive{archive_stream};
        p_jst_original.save(archive);
    }

    libjst::journal_sequence_tree_partitioned p_jst_copy{std::addressof(jst)};
    // deserialise
    {
        cereal::BinaryInputArchive archive{archive_stream};
        p_jst_copy.load(archive);
    }

    EXPECT_EQ(p_jst_original.bin_count(), p_jst_copy.bin_count());

    for (uint32_t index = 0; index < p_jst_original.bin_count(); ++index)
    {
        auto context_enumerator1 = p_jst_original.context_enumerator(libjst::context_size{GetParam().context_size},
                                                                     libjst::bin_index{index});
        auto context_enumerator2 = p_jst_copy.context_enumerator(libjst::context_size{GetParam().context_size},
                                                                 libjst::bin_index{index});

        auto context_it1 = context_enumerator1.begin();
        auto context_it2 = context_enumerator2.begin();
        for (; context_it1 != context_enumerator1.end(); ++context_it1, ++context_it2)
        {
            EXPECT_RANGE_EQ(*context_it1, *context_it2);
            EXPECT_RANGE_EQ(context_it1.positions(), context_it2.positions());
        }
    }
}

// ----------------------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(substitution_with_one_bin, partitioned_traversal_test, testing::Values(
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
    .reference{"aaaabbb"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{4u},
    .bin_count{1u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_with_second_bin_empty, partitioned_traversal_test, testing::Values(
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
    //              c
    // bin0: "aaaa|bbb"   aaaa, aaab, aabc, abcb, aabb, abbb
    // bin1: "    |bbb"; // empty.
    .reference{"aaaabbb"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_with_two_bins, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //               c
    //          0123456
    // bin0:   "aaaa"   : aa, aa, aa, ab
    // bin1:   "    bbb": bc, cb, bb, bb
    .reference{"aaaabbb"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{2u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_on_boundary, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //         "     c
    // bin0:   "aaaa|   ": aaa, aaa, aac, acb, aaa, aaa
    // bin1:   "    |bbb": cbb, bbb
    .reference{"aaaabbb"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{4u, libjst::test::substitution_t{"c"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_ending_in_boundary, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //         "   b
    // bin0:   "aaaa"   : aaa, aab, aba, baa, aaa, aaa, aaa
    // bin1:   "    aaa": aaa
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{3u, libjst::test::substitution_t{"b"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(insertion_on_boundary, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //         "     bb
    // bin0:   "aaaa|  "   : aaa, aaa, aab, abb, aaa, aaa
    // bin1:   "    |  aaa": bba, baa, aaa
    .reference{"aaaaaaa"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{4u, libjst::test::insertion_t{"bb"s}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(deletion_over_boundary, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //         "
    // bin0:   "aaa-   ": aaa, aa--b, a--bb, aaa, aab, abb
    // bin1:   "    -bb": bbb [0: 4, 3: 4]
    //          aaaabb|
    //             --
    //          0123456
    .reference{"aaaabbb"s},
    .sequence_count{4u},
    .events
    {
        libjst::test::shared_event_t{3u, libjst::test::deletion_t{2}, libjst::test::coverage_t{0, 1, 1, 0}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(insertion_before_bin_boundary, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s},
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{4u, libjst::test::insertion_t{"iii"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(insertion_at_end_of_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s},
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{8u, libjst::test::insertion_t{"iii"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_of_entire_last_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s}, // |r| = 10
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::substitution_t{"ccccc"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(substitution_of_entire_first_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s}, // |r| = 10
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{0u, libjst::test::substitution_t{"ccccc"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(deletion_of_entire_last_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s}, // |r| = 10
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{5u, libjst::test::deletion_t{5}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(deletion_of_entire_first_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    .reference{"aaaaabbbbb"s}, // |r| = 10
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{0u, libjst::test::deletion_t{5}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{3u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(deletion_overlaps_bin_and_substitution_in_second_bin, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //        0123456789
    //        aaaaabbbbb
    // s0:    aaaa_____b
    // s1:    aaaaabbbcc

    .reference{"aaaaabbbbb"s},
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{4u, libjst::test::deletion_t{5}, libjst::test::coverage_t{1, 0}},
        libjst::test::shared_event_t{8u, libjst::test::substitution_t{"cc"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(insertion_with_substitution_on_last_position, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //        01234567890123
    //        012345    6789
    //        aaaaab----cccc
    // s0:    aaaaabkkkkcccc
    // s1:    aaaaab----ccrr
    .reference{"aaaaabcccc"s},
    .sequence_count{2u},
    .events
    {
        libjst::test::shared_event_t{6u, libjst::test::insertion_t{"kkkk"s}, libjst::test::coverage_t{1, 0}},
        libjst::test::shared_event_t{8u, libjst::test::substitution_t{"rr"s}, libjst::test::coverage_t{0, 1}}
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(complex_tree_with_two_bins, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //        0123    45    6789
    //        xaaa----bb----cccy
    // s0: f--x---ii--qq----qqqylll
    // s1: gg-ppppii--bbkkkkcccylll
    // s2: hhhx---jjjjqq----qqqy
    // s3:    ppppjjjjbbkkkkcccy
    // s4:    xaaaii--__----___y
    // s5:    xaaaii--bb----cccy
    // s6:    xaaajjjj__----___ylll
    // s7:    xaaajjjjbb----ccrrlll

    .reference{"xaaabbcccy"s},
    .sequence_count{8u},
    .events
    {
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"f"s},
                                     libjst::test::coverage_t{1, 0, 0, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"gg"s},
                                     libjst::test::coverage_t{0, 1, 0, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"hhh"s},
                                     libjst::test::coverage_t{0, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::substitution_t{"pppp"s},
                                     libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
        libjst::test::shared_event_t{1u,
                                     libjst::test::deletion_t{3},
                                     libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::insertion_t{"ii"s},
                                     libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::insertion_t{"jjjj"s},
                                     libjst::test::coverage_t{0, 0, 1, 1, 0, 0, 1, 1}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::substitution_t{"qqqqq"s},
                                     libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::deletion_t{5},
                                     libjst::test::coverage_t{0, 0, 0, 0, 1, 0, 1, 0}},
        libjst::test::shared_event_t{6u,
                                     libjst::test::insertion_t{"kkkk"s},
                                     libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
        libjst::test::shared_event_t{8u,
                                     libjst::test::substitution_t{"rr"s},
                                     libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 1}},
        libjst::test::shared_event_t{10u,
                                     libjst::test::insertion_t{"lll"s},
                                     libjst::test::coverage_t{1, 1, 0, 0, 0, 1, 0, 1}},
    },
    .context_size{4u},
    .bin_count{2u}
}));

INSTANTIATE_TEST_SUITE_P(complex_tree_with_three_bins, partitioned_traversal_test, testing::Values(
libjst::test::traversal_fixture
{
    //        0123    45    6789
    //        xaaa----bb----cccy
    // s0: f--x---ii--qq----qqqylll
    // s1: gg-ppppii--bbkkkkcccylll
    // s2: hhhx---jjjjqq----qqqy
    // s3:    ppppjjjjbbkkkkcccy
    // s4:    xaaaii--__----___y
    // s5:    xaaaii--bb----cccy
    // s6:    xaaajjjj__----___ylll
    // s7:    xaaajjjjbb----ccrrlll

    .reference{"xaaabbcccy"s},
    .sequence_count{8u},
    .events
    {
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"f"s},
                                     libjst::test::coverage_t{1, 0, 0, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"gg"s},
                                     libjst::test::coverage_t{0, 1, 0, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::insertion_t{"hhh"s},
                                     libjst::test::coverage_t{0, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{0u,
                                     libjst::test::substitution_t{"pppp"s},
                                     libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
        libjst::test::shared_event_t{1u,
                                     libjst::test::deletion_t{3},
                                     libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::insertion_t{"ii"s},
                                     libjst::test::coverage_t{1, 1, 0, 0, 1, 1, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::insertion_t{"jjjj"s},
                                     libjst::test::coverage_t{0, 0, 1, 1, 0, 0, 1, 1}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::substitution_t{"qqqqq"s},
                                     libjst::test::coverage_t{1, 0, 1, 0, 0, 0, 0, 0}},
        libjst::test::shared_event_t{4u,
                                     libjst::test::deletion_t{5},
                                     libjst::test::coverage_t{0, 0, 0, 0, 1, 0, 1, 0}},
        libjst::test::shared_event_t{6u,
                                     libjst::test::insertion_t{"kkkk"s},
                                     libjst::test::coverage_t{0, 1, 0, 1, 0, 0, 0, 0}},
        libjst::test::shared_event_t{8u,
                                     libjst::test::substitution_t{"rr"s},
                                     libjst::test::coverage_t{0, 0, 0, 0, 0, 0, 0, 1}},
        libjst::test::shared_event_t{10u,
                                     libjst::test::insertion_t{"lll"s},
                                     libjst::test::coverage_t{1, 1, 0, 0, 0, 1, 0, 1}},
    },
    .context_size{2u},
    .bin_count{3u}
}));
