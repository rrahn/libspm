// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <libjst/sequence/journal_entry.hpp>

// Scenario: Initialization of journal entry
// Given a default initialized journal entry
// When accessing the begin position
// Then the begin position is 0
// And the segment is empty
// And the end position is 0
TEST(journal_entry_test, default_constructor)
{
    using journal_entry_t = libjst::journal_entry<size_t, std::vector<char>::iterator>;
    EXPECT_TRUE((std::default_initializable<journal_entry_t>));
    journal_entry_t entry{};
    EXPECT_EQ(entry.begin_position(), 0u);
    EXPECT_EQ(entry.end_position(), 0u);
    EXPECT_TRUE(entry.segment().empty());
}

// Given a journal entry initialized with a position and a sequence segment
// When accessing the begin position
// Then the begin position is returned
TEST(journal_entry_test, begin_position)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry{42u, sequence};
    EXPECT_EQ(entry.begin_position(), 42u);
}

// Given a journal entry initialized with a position and a sequence segment
// When accessing the end position
// Then the begin position + the size of the segment is returned
TEST(journal_entry_test, end_position)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry{42u, sequence};
    EXPECT_EQ(entry.end_position(), 46u);
}

// When accessing the segment
// Then the segment is returned
TEST(journal_entry_test, segment)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry{42u, sequence};
    EXPECT_EQ(entry.segment().data(), sequence.data());
}

// Scenario: Compare two journal entries
// Given two initialized journal entries
// When both journal entries have the same begin position and the same segment
// Then the journal entries are equivalent
TEST(journal_entry_test, equality)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry1{42u, sequence};
    libjst::journal_entry entry2{42u, sequence};
    EXPECT_EQ(entry1 <=> entry2, std::weak_ordering::equivalent);
    EXPECT_EQ(entry1, entry2);
}

// Given two initialized journal entries
// When both journal entries have different begin positions but the same segment
// Then the journal entries are not equivalent
TEST(journal_entry_test, inequality_begin_position)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry1{42u, sequence};
    libjst::journal_entry entry2{43u, sequence};
    EXPECT_NE(entry1, entry2);
}

// Given two initialized journal entries
// When both journal entries hold different segments but the same begin position
// Then the journal entries are equivalent
TEST(journal_entry_test, inequality_segment)
{
    std::vector sequence1{'A', 'C', 'G', 'T'};
    std::vector sequence2{'A', 'C', 'G', 'T', 'A'};
    libjst::journal_entry entry1{42u, sequence1};
    libjst::journal_entry entry2{42u, sequence2};
    EXPECT_EQ(entry1, entry2);
}

// Given two initialized journal entries
// When both journal entries hold different segments and different begin positions
// Then the journal entries are not equivalent
TEST(journal_entry_test, inequality_segment_and_begin_position)
{
    std::vector sequence1{'A', 'C', 'G', 'T'};
    std::vector sequence2{'A', 'C', 'G', 'T', 'A'};
    libjst::journal_entry entry1{42u, sequence1};
    libjst::journal_entry entry2{43u, sequence2};
    EXPECT_NE(entry1, entry2);
}

// Scenario: Compare ordering of two journal entries
// Given two initialized journal entries
// When the begin position of the first journal entry is less than the begin position of the second journal entry
// Then the first journal entry is less than the second journal entry
TEST(journal_entry_test, less_than)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry1{42u, sequence};
    libjst::journal_entry entry2{43u, sequence};
    EXPECT_EQ(entry1 <=> entry2, std::weak_ordering::less);
    EXPECT_LT(entry1, entry2);
}

// Given two initialized journal entries
// When the begin position of the first journal entry is greater than the begin position of the second journal entry
// Then the first journal entry is greater than the second journal entry
TEST(journal_entry_test, greater_than)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journal_entry entry1{43u, sequence};
    libjst::journal_entry entry2{42u, sequence};
    EXPECT_EQ(entry1 <=> entry2, std::weak_ordering::greater);
    EXPECT_GT(entry1, entry2);
}
