// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/journal.hpp>

struct journal_test : public ::testing::Test
{
    using key_type = uint32_t;
    using mapped_type = libjst::detail::subrange_t<std::string &>;
    using journal_type = libjst::journal<key_type, std::string &>;

    std::string sequence{"aaaaccccggggtttt"};
};

// TODO: check iterator concept.

TEST_F(journal_test, construction)
{
    journal_type journal{};
    EXPECT_TRUE(journal.empty());

    journal = journal_type{sequence};
    EXPECT_FALSE(journal.empty());
}

TEST_F(journal_test, sequence)
{
    journal_type journal{sequence};
    EXPECT_RANGE_EQ(journal.sequence(), sequence);

    // TODO: Check const iterator
}

TEST_F(journal_test, record_insertion)
{
    using namespace std::literals;

    std::string segment{"uu"s};

    { // insert in middel
        journal_type journal{sequence};

        auto [pos, seq] = *(journal.record_insertion(8u, segment));
        EXPECT_EQ(pos, 8u);
        EXPECT_RANGE_EQ(seq, segment);
        auto journal_range = journal.sequence();

        EXPECT_EQ(journal_range.size(), sequence.size() + segment.size());
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal_range, expected.insert(8, segment));
    }

    // insert in empty
    {
        // expect throw or false as return type?
        journal_type journal{};
        // TODO: return pair with false and end iterator?
        // auto [pos, seq] = *journal.record_insertion(8, segment);

        // EXPECT_FALSE(insert_result);
        // EXPECT_RANGE_EQ(seq, segment);
        // EXPECT_RANGE_EQ(journal.sequence(), std::string{}); TODO: check me!

        auto [pos, seq] = *journal.record_insertion(0, segment);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, segment);
        EXPECT_RANGE_EQ(journal.sequence(), segment);
    }

    // insert at end
    {
        journal_type journal{sequence};
        // auto [pos, seq] = *journal.record_insertion(sequence.size() + 1, segment);

        // EXPECT_FALSE(insert_result);
        // EXPECT_RANGE_EQ(seq, segment.size());
        // EXPECT_RANGE_EQ(journal.sequence(), sequence);

        auto [pos, seq] = *journal.record_insertion(sequence.size(), segment);
        EXPECT_EQ(pos, sequence.size());
        EXPECT_RANGE_EQ(seq, segment);
        EXPECT_RANGE_EQ(journal.sequence(), sequence + segment);
    }

    // insert at beginning
    {
        journal_type journal{sequence};
        auto [pos, seq] = *journal.record_insertion(0, segment);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, segment);
        EXPECT_RANGE_EQ(journal.sequence(), segment + sequence);
    }

    // insert on same position twice.
    {
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_insertion(8, segment);
        EXPECT_EQ(pos, 8u);
        EXPECT_RANGE_EQ(seq, segment);

        auto [pos2, seq2] = *journal.record_insertion(8, segment);
        EXPECT_EQ(pos2, 8u);
        EXPECT_RANGE_EQ(seq2, segment);
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.insert(8, segment).insert(8, segment));
    }
}

TEST_F(journal_test, record_insertion_in_empty_journal_sequence)
{
    using namespace std::literals;

    std::string empty_ref{""s};
    std::string single_insertion{"i"s};

    journal_type journal{empty_ref};
    EXPECT_FALSE(journal.empty());
    EXPECT_RANGE_EQ(journal.sequence(), empty_ref);

    auto [pos, seq] = *journal.record_insertion(0, single_insertion);

    EXPECT_EQ(pos, 0);
    EXPECT_RANGE_EQ(seq, single_insertion);
    EXPECT_FALSE(journal.empty());
    EXPECT_RANGE_EQ(journal.sequence(), single_insertion);
}

TEST_F(journal_test, record_deletion)
{
    using namespace std::literals;

    //-------------------------------------------------------------------------
    // invalid erase

    // { // erase from empty journal decorator
    //     journal_type journal{};

    //     auto [pos, seq] = *journal.record_deletion(0, 10);

    //     EXPECT_EQ(pos, 0);
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_TRUE(journal.empty());
    // }

    // { // erase inverted range
    //     journal_type journal{sequence};

    //     auto [pos, seq] = *journal.record_deletion(5, 4);

    //     EXPECT_EQ(pos, );
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    // { // erase region that goes over the end of the journal decorator
    //     journal_type journal{sequence};

    //     auto [pos, seq] = *journal.record_deletion(5, 17);

    //     EXPECT_EQ(pos, );
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    // { // erase region that goes over the end of the journal decorator
    //     journal_type journal{sequence};

    //     auto [pos, seq] = *journal.record_deletion(16, 17);

    //     EXPECT_EQ(pos, );
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    // { // erase empty segment
    //     journal_type journal{sequence};

    //     auto [pos, seq] = *journal.record_deletion(5, 5);

    //     EXPECT_EQ(pos, );
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    //-------------------------------------------------------------------------
    // erase from journal decorator with single entry

    { // erase some internal segment
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_deletion(4, 4);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, sequence.substr(8));
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(4, 4));
    }

    { // erase single element within entry
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_deletion(7, 1);

        EXPECT_EQ(pos, 7u);
        EXPECT_RANGE_EQ(seq, sequence.substr(8));
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(7, 1));
    }

    { // erase entire entry
        journal_type journal{sequence};

        auto it = journal.record_deletion(0, 16);
        EXPECT_EQ(it, std::ranges::end(journal));
        EXPECT_RANGE_EQ(journal.sequence(), ""s);
    }

    { // erase suffix of entry
        journal_type journal{sequence};

        auto it = journal.record_deletion(5, 11);
        EXPECT_EQ(it, std::ranges::end(journal));
        EXPECT_EQ(journal.size(), 1u);
        EXPECT_RANGE_EQ(journal.sequence(), sequence.substr(0, 5));
    }

    // erase prefix of entry
    {
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_deletion(0, 5);

        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, sequence.substr(5));
        EXPECT_EQ(journal.size(), 1u);
        EXPECT_RANGE_EQ(journal.sequence(), sequence.substr(5));
    }

    //-------------------------------------------------------------------------
    // erase from journal decorator with multiple entries
    // erase between two adjacent entries

    journal_type journal_base{sequence};

    EXPECT_NE((journal_base.record_deletion(12, 1)), std::ranges::end(journal_base));
    EXPECT_NE((journal_base.record_deletion(8, 1)), std::ranges::end(journal_base));
    EXPECT_NE((journal_base.record_deletion(4, 1)), std::ranges::end(journal_base));
    EXPECT_NE((journal_base.record_deletion(0, 1)), std::ranges::end(journal_base));

    std::string expected_base{sequence};
    EXPECT_RANGE_EQ(journal_base.sequence(), expected_base.erase(12, 1).erase(8, 1).erase(4, 1).erase(0, 1));

    { // erase begin from first to last
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(3, 6);
        EXPECT_EQ(pos, 3u);
        EXPECT_RANGE_EQ(seq, "ttt"s);

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(3, 6));
    }

    { // from last of first to first of last // two elements are delted
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(5, 2);
        EXPECT_EQ(pos, 5u);
        EXPECT_RANGE_EQ(seq, "gg"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(5, 2));
    }

    { // from middle to middle - some range is deleted
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(4, 4);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, "g"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(4, 4));
    }

    { // from middle to end of second - entire second entry is deleted
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(4, 5);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, "ttt"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(4, 5));
    }

    { // from begin to middle of second - entire first entry is deleted
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(3, 5);
        EXPECT_EQ(pos, 3u);
        EXPECT_RANGE_EQ(seq, "g"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(3, 5));
    }

    //-------------------------------------------------------------------------
    // erase from journal decorator with multiple entries
    // erase between two distant entries

    { // erase begin from first to last
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(0, 9);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, "ttt"s);
        EXPECT_FALSE(journal.empty());
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(0, 9));
    }

    { // erase everything
        journal_type journal{journal_base};

        auto it = journal.record_deletion(0, 12);
        EXPECT_EQ(it, std::ranges::end(journal));
        EXPECT_RANGE_EQ(journal.sequence(), ""s);
    }

    { // from last of first to first of last - two elements are delted and all entries in between
        journal_type journal{journal_base};
        auto [pos, seq] = *journal.record_deletion(2, 8);
        EXPECT_EQ(pos, 2u);
        EXPECT_RANGE_EQ(seq, "tt"s);
        EXPECT_FALSE(journal.empty());
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(2, 8));
    }

    { // from last of first to first of last - two elements are delted and all entries in between
        journal_type journal{journal_base};
        auto [pos, seq] = *journal.record_deletion(1, 9);
        EXPECT_EQ(pos, 1u);
        EXPECT_RANGE_EQ(seq, "tt"s);
        EXPECT_FALSE(journal.empty());
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(1, 9));
    }

    { // from middle to middle - some elements are delted and all entries in between
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(1, 10);
        EXPECT_EQ(pos, 1u);
        EXPECT_RANGE_EQ(seq, "t"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(1, 10));
    }

    { // from middle to end - some elements are delted and all entries in between
        journal_type journal{journal_base};

        auto it = journal.record_deletion(1, 11);
        EXPECT_EQ(it, std::ranges::end(journal));
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(1, 11));
    }

    { // from begin to middle - some elements are delted and all entries in between
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_deletion(0, 11);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, "t"s);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.erase(0, 11));
    }
}

TEST_F(journal_test, record_substitution)
{
    using namespace std::literals;

    std::string segment{"uu"};

    //-------------------------------------------------------------------------
    // replace invalid

    // { // replace in empty journal decorator
    //     journal_type journal{};

    //     EXPECT_FALSE(journal.record_substitution(0, segment));
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_TRUE(journal.empty());

    //     EXPECT_FALSE(journal.record_substitution(10, segment));
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_TRUE(journal.empty());
    // }

    // { // position is after end
    //     journal_type journal{sequence};

    //     EXPECT_FALSE(journal.record_substitution(17, segment));
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    // { // last position is after end
    //     journal_type journal{sequence};

    //     EXPECT_FALSE(journal.record_substitution(15, segment));
    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    // { // segment is empty
    //     journal_type journal{sequence};

    //     EXPECT_FALSE(journal.record_substitution(10, mapped_type}));    //     EXPECT_RANGE_EQ(seq, segment);
    //     EXPECT_FALSE(journal.empty());
    // }

    //-------------------------------------------------------------------------
    // replace within single entry

    { // replace middle of entry
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_substitution(4, segment);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, segment);
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(4, segment.size(), segment));
    }

    { // replace entire entry
        journal_type journal{sequence};

        std::string replace_all(sequence.size(), 'u');
        auto [pos, seq] = *journal.record_substitution(0, replace_all);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, replace_all);
        EXPECT_RANGE_EQ(journal.sequence(), replace_all);
    }

    { // replace prefix of entry
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_substitution(0, segment);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, segment);
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(0, segment.size(), segment));
    }

    { // replace suffix of entry
        journal_type journal{sequence};

        auto [pos, seq] = *journal.record_substitution(14, segment);
        EXPECT_EQ(pos, 14u);
        EXPECT_RANGE_EQ(seq, segment);
        std::string expected{sequence};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(14, segment.size(), segment));
    }

    //-------------------------------------------------------------------------
    // replace with multiple entries
    // replace between adjacent entries

    journal_type journal_base{sequence};
    // create: aa|uu|cc|uu|gg|uu|tt|uu
    journal_base.record_substitution(2, segment);
    journal_base.record_substitution(6, segment);
    journal_base.record_substitution(10, segment);
    journal_base.record_substitution(14, segment);

    std::string expected_base{sequence};
    expected_base.replace(2, segment.size(), segment).
                  replace(6, segment.size(), segment).
                  replace(10, segment.size(), segment).
                  replace(14, segment.size(), segment);

    {
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(5, segment);
        EXPECT_EQ(pos, 5u);
        EXPECT_RANGE_EQ(seq, segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(5, segment.size(), segment));
    }

    {
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(8, segment);
        EXPECT_EQ(pos, 8u);
        EXPECT_RANGE_EQ(seq, segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(8, segment.size(), segment));
    }

    { // replace from begin of node to end of next node.
        std::string new_segment{"xxxx"};
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(4, new_segment);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(4, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"xxx"};
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(4, new_segment);
        EXPECT_EQ(pos, 4u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(4, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"xxx"};
        journal_type journal{journal_base};
        auto [pos, seq] = *journal.record_substitution(5, new_segment);
        EXPECT_EQ(pos, 5u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(5, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"x"};
        journal_type journal{journal_base};
        auto [pos, seq] = *journal.record_substitution(5, new_segment);
        EXPECT_EQ(pos, 5u);
        EXPECT_RANGE_EQ(seq, new_segment);
        auto [pos2, seq2] = *journal.record_substitution(4, new_segment);
        EXPECT_EQ(pos2, 4u);
        EXPECT_RANGE_EQ(seq2, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(5, new_segment.size(), new_segment).
                                     replace(4, new_segment.size(), new_segment));
    }

    //-------------------------------------------------------------------------
    // replace with multiple entries
    // replace between distant entries

    { // entire range
        std::string new_segment(sequence.size(), 'y');
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(0, new_segment);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, new_segment);
        EXPECT_RANGE_EQ(journal.sequence(), new_segment);
    }

    std::string new_segment(sequence.size() - 5, 'y');
    { // replace prefix
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(0, new_segment);
        EXPECT_EQ(pos, 0u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(0, new_segment.size(), new_segment));
    }

    { // replace suffix
        journal_type journal{journal_base};
        size_t const replace_position = journal.sequence().size() - new_segment.size();
        auto [pos, seq] = *journal.record_substitution(replace_position, new_segment);
        EXPECT_EQ(pos, replace_position);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(replace_position, new_segment.size(), new_segment));
    }

    { // replace middle from begin
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(3, new_segment);
        EXPECT_EQ(pos, 3u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(3, new_segment.size(), new_segment));
    }

    {
        journal_type journal{journal_base};

        auto [pos, seq] = *journal.record_substitution(2, new_segment);
        EXPECT_EQ(pos, 2u);
        EXPECT_RANGE_EQ(seq, new_segment);
        std::string expected{expected_base};
        EXPECT_RANGE_EQ(journal.sequence(), expected.replace(2, new_segment.size(), new_segment));
    }
}
