// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <libjst/sequence/journaled_sequence.hpp>

// Scenario: Initialization of journaled_sequence
// Given a default initialized journaled_sequence
// Then the journaled_sequence is empty
// And its size is 0
TEST(journaled_sequence_test, default_constructor)
{
    using journaled_sequence_t = libjst::journaled_sequence<std::vector<char> &, size_t>;
    EXPECT_TRUE((std::default_initializable<journaled_sequence_t>));
    journaled_sequence_t journaled_sequence{};
    EXPECT_TRUE(journaled_sequence.empty());
    EXPECT_EQ(journaled_sequence.size(), 0u);
}

// Given a journaled_sequence referencing memory from std::vector
// Then the journaled_sequence satisfies the constraints for a random access range
TEST(journaled_sequence_test, random_access_range)
{
    using journaled_sequence_t = libjst::journaled_sequence<std::vector<char> &, size_t>;
    EXPECT_TRUE((std::ranges::random_access_range<journaled_sequence_t>));
}

// Given a default initialized journaled_sequence
// When inserting a non-empty sequence segment
// Then the journaled_sequence is not empty and its range spells the inserted sequence
TEST(journaled_sequence_test, insert_empty)
{
    using journaled_sequence_t = libjst::journaled_sequence<std::vector<char> &, size_t>;
    journaled_sequence_t journaled_sequence{};
    std::vector sequence{'A', 'C', 'G', 'T'};
    journaled_sequence.insert(journaled_sequence.begin(), sequence);
    EXPECT_FALSE(journaled_sequence.empty());
    EXPECT_EQ(journaled_sequence.size(), sequence.size());
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, sequence));
}

// Given an initialized journaled_sequence
// When inserting a non-empty sequence segment into the middle of the journaled_sequence
// When inserting a second sequence segment before the first inserted sequence
// When inserting a third sequence segment after the first inserted sequence
// Then the prefix and suffix of the journaled_sequence are not changed and the inserted sequences are inserted at the correct positions
TEST(journaled_sequence_test, insert_middle)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    std::vector sequence1{'T', 'G', 'C', 'A'};
    auto it = journaled_sequence.insert(journaled_sequence.begin() + 2, sequence1);
    std::vector sequence2{'G', 'G'};
    it = journaled_sequence.insert(it, sequence2);
    std::vector sequence3{'C', 'C'};
    journaled_sequence.insert(it + 6, sequence3);
    std::vector expected_sequence{'A', 'C', 'G', 'G', 'T', 'G', 'C', 'A', 'C', 'C', 'G', 'T'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence));
}

// Given an initialized journaled_sequence
// When erasing an empty segment from it
// Then the journaled_sequence is not changed
// And the returned iterator is the same as the end of the erased segment
TEST(journaled_sequence_test, erase_empty_segment)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    auto erase_position = journaled_sequence.begin() + 2;
    auto it = journaled_sequence.erase(erase_position, erase_position);
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, sequence));
    EXPECT_EQ(it, erase_position);
}

// Given an initialized journaled_sequence
// When erasing a non-empty segment from the middle of it
// When erasing a second single element before the first erased segment
// When erasing a third single element after the first erased segment
// When erasing a fourth segment that spans the whole journaled_sequence
// Then first the prefix and suffix of the journaled_sequence are not changed and the erased segment is erased
// And the single elements are erased
// And finally the whole journaled_sequence is erased
TEST(journaled_sequence_test, erase_nested_segments)
{
    std::vector sequence{'A', 'A', 'C', 'C', 'G', 'G', 'T', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    auto erase_position = journaled_sequence.begin() + 3;
    auto it = journaled_sequence.erase(erase_position, erase_position + 2);
    std::vector expected_sequence1{'A', 'A', 'C', 'G', 'T', 'T'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence1));
    EXPECT_EQ(journaled_sequence.size(), expected_sequence1.size());

    it = journaled_sequence.erase(it - 1);
    std::vector expected_sequence2{'A', 'A', 'G', 'T', 'T'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence2));
    EXPECT_EQ(journaled_sequence.size(), expected_sequence2.size());

    it = journaled_sequence.erase(it);
    std::vector expected_sequence3{'A', 'A', 'T', 'T'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence3));
    EXPECT_EQ(journaled_sequence.size(), expected_sequence3.size());

    it = journaled_sequence.erase(journaled_sequence.begin(), journaled_sequence.end());
    EXPECT_TRUE(journaled_sequence.empty());
    EXPECT_EQ(it, journaled_sequence.end());
}

// Given an initialized journaled_sequence
// When replacing a non-empty segment with a non-empty sequence segment
// Then the journaled_sequence spells the new sequence
TEST(journaled_sequence_test, replace_segment)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    std::vector sequence1{'T', 'G', 'C', 'A'};
    auto it = journaled_sequence.replace(journaled_sequence.begin() + 2, journaled_sequence.begin() + 4, sequence1);
    std::vector expected_sequence{'A', 'C', 'T', 'G', 'C', 'A'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence));
    EXPECT_EQ(it, journaled_sequence.begin() + 2);
    EXPECT_EQ(journaled_sequence.size(), expected_sequence.size());
}

// Given an initialized journaled sequence
// When replacing a non-empty segment with an empty sequence segment
// Then the journaled_sequence is not changed
TEST(journaled_sequence_test, replace_segment_with_empty)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    std::vector<char> sequence1{};
    auto it = journaled_sequence.replace(journaled_sequence.begin() + 2, journaled_sequence.begin() + 4, sequence1);
    std::vector expected_sequence{'A', 'C'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence));
    EXPECT_EQ(it, journaled_sequence.begin() + 2);
    EXPECT_EQ(journaled_sequence.size(), expected_sequence.size());
}

// Given an initialized journaled_sequence
// When replacing an empty segment with a non-empty sequence segment
// Then the journaled_sequence is modified as if the segment was inserted behind the end of the erased.
TEST(journaled_sequence_test, replace_empty_segment)
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    std::vector sequence1{'T', 'G', 'C', 'A'};
    auto it = journaled_sequence.replace(journaled_sequence.begin() + 2, journaled_sequence.begin() + 2, sequence1);
    std::vector expected_sequence{'A', 'C', 'T', 'G', 'C', 'A', 'G', 'T'};
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence));
    EXPECT_EQ(it, journaled_sequence.begin() + 2);
    EXPECT_EQ(journaled_sequence.size(), expected_sequence.size());
}

// Given an initialized journaled_sequence
// When the journaled_sequence is modified using the insert, erase, and replace functions with const iterators
// Then the resulting journaled_sequence spells the expected sequence with modifications
// And the original sequence is not changed
TEST(journaled_sequence_test, modify_using_const_iterators)
{
    using namespace std::literals;

    std::vector sequence{'A', 'C', 'G', 'T'};
    libjst::journaled_sequence journaled_sequence{sequence};
    std::vector sequence1{'T', 'G', 'C', 'A'};
    auto it = journaled_sequence.insert(std::as_const(journaled_sequence).begin() + 2, sequence1);

    auto erase_position = std::as_const(journaled_sequence).begin() + 3;
    it = journaled_sequence.erase(erase_position, erase_position + 2);
    std::vector sequence2{'G', 'G'};
    it = journaled_sequence.replace(std::as_const(journaled_sequence).begin() + 2,
                                    std::as_const(journaled_sequence).begin() + 4, sequence2);

    std::vector expected_sequence{'A', 'C', 'G', 'G', 'G', 'T'};
    EXPECT_EQ(journaled_sequence.size(), expected_sequence.size());
    EXPECT_TRUE(std::ranges::equal(journaled_sequence, expected_sequence));
    EXPECT_TRUE(std::ranges::equal(sequence, "ACGT"s));
}
