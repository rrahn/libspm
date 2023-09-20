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
