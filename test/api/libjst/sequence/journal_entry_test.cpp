// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include <libjst/sequence/journal_entry.hpp>

SCENARIO("Journal entry can be default initialized", "[sequence][journal_entry]")
{
    GIVEN("A default initialized journal entry")
    {
        using journal_entry_t = libjst::journal_entry<size_t, std::vector<char>::iterator>;
        journal_entry_t entry{};
        THEN("the begin position is 0")
        {
            REQUIRE(entry.begin_position() == 0u);
        }
        AND_THEN("the segment is empty")
        {
            REQUIRE(entry.segment().empty());
        }
        AND_THEN("the end position is 0")
        {
            REQUIRE(entry.end_position() == 0u);
        }
    }
}

SCENARIO("Journal entry can be initialized with a position and a sequence", "[sequence],[journal_entry]")
{
    GIVEN("A journal entry initialized with a position and a sequence segment")
    {
        std::vector sequence{'A', 'C', 'G', 'T'};
        libjst::journal_entry entry{42u, sequence};
        THEN("the begin position is the given position")
        {
            REQUIRE(entry.begin_position() == 42u);
        }
        AND_THEN("the segment is the given sequence")
        {
            REQUIRE(entry.segment().data() == sequence.data());
        }
        AND_THEN("the end position is the given position + the size of the sequence")
        {
            REQUIRE(entry.end_position() == 46u);
        }
    }
}

SCENARIO("Two journal entries can be compared", "[sequence],[journal_entry]")
{
    std::vector sequence{'A', 'C', 'G', 'T'};
    GIVEN("Two journal entries with the same begin position and the same segment")
    {
        libjst::journal_entry entry1{42u, sequence};
        libjst::journal_entry entry2{42u, sequence};
        THEN("the journal entries are equivalent")
        {
            REQUIRE(entry1 <=> entry2 == std::weak_ordering::equivalent);
            REQUIRE(entry1 == entry2);
        }
    }
    AND_GIVEN("Two journal entries with the same begin position but different segments")
    {
        std::vector sequence2{'A', 'C', 'G', 'T'};
        libjst::journal_entry entry1{42u, sequence};
        libjst::journal_entry entry2{42u, sequence2};
        THEN("the journal entries are not equivalent")
        {
            REQUIRE(entry1 != entry2);
        }
    }
    AND_GIVEN("Two journal entries with different begin positions but the same segment")
    {
        libjst::journal_entry entry1{42u, sequence};
        libjst::journal_entry entry2{43u, sequence};
        THEN("the journal entries are not equivalent")
        {
            REQUIRE(entry1 != entry2);
        }
    }
    AND_GIVEN("A journal entry with lower begin position than the second journal entry")
    {
        libjst::journal_entry entry1{42u, sequence};
        libjst::journal_entry entry2{43u, sequence};
        THEN("the first journal entry is less than the second journal entry")
        {
            REQUIRE(entry1 <=> entry2 == std::weak_ordering::less);
            REQUIRE(entry1 < entry2);
        }
        AND_THEN("the second journal entry is greater than the first journal entry")
        {
            REQUIRE(entry2 <=> entry1 == std::weak_ordering::greater);
            REQUIRE(entry2 > entry1);
        }
    }
}

SCENARIO("A journal entry can be split")
{
    GIVEN("A journal entry with begin position 42 and segment 'ACGT'")
    {
        std::vector sequence{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
        auto segment_begin = sequence.begin() + 4;
        auto segment_end = sequence.begin() + 8;
        libjst::journal_entry entry{42u, std::ranges::subrange{segment_begin, segment_end}};
        WHEN("the journal entry is split at position 44")
        {
            auto split_it = entry.segment().begin() + 2;
            auto [entry1, entry2] = entry.split_at(split_it);
            THEN("the first journal entry covers the interval [42, 44)")
            {
                REQUIRE(entry1.begin_position() == 42u);
                REQUIRE(entry1.segment().begin() == segment_begin);
                REQUIRE(entry1.end_position() == 44u);
            }
            AND_THEN("the second journal entry covers the interval [44, 46)")
            {
                REQUIRE(entry2.begin_position() == 44u);
                REQUIRE(entry2.segment().begin() == segment_begin + 2);
                REQUIRE(entry2.end_position() == 46u);
            }
        }
        WHEN("the journal entry is split before the begin of the segment")
        {
            auto split_it = entry.segment().begin() - 1;
            auto [entry1, entry2] = entry.split_at(split_it);
            THEN("the first journal entry covers the interval [42, 42)")
            {
                REQUIRE(entry1.begin_position() == 42u);
                REQUIRE(entry1.segment().begin() == segment_begin);
                REQUIRE(entry1.end_position() == 42u);
            }
            AND_THEN("the second journal entry covers the interval [42, 46)")
            {
                REQUIRE(entry2.begin_position() == 42u);
                REQUIRE(entry2.segment().begin() == segment_begin);
                REQUIRE(entry2.end_position() == 46u);
            }
        }
        WHEN("the journal entry is split after the end of the segment")
        {
            auto split_it = entry.segment().end() + 1;
            auto [entry1, entry2] = entry.split_at(split_it);
            THEN("the first journal entry covers the interval [42, 46)")
            {
                REQUIRE(entry1.begin_position() == 42u);
                REQUIRE(entry1.segment().begin() == segment_begin);
                REQUIRE(entry1.end_position() == 46u);
            }
            AND_THEN("the second journal entry covers the interval [46, 46)")
            {
                REQUIRE(entry2.begin_position() == 46u);
                REQUIRE(entry2.segment().begin() == segment_end);
                REQUIRE(entry2.end_position() == 46u);
            }
        }
    }
}
