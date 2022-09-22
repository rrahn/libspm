// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <ranges>

#include <seqan3/test/expect_range_eq.hpp>

#include <libio/file/pivot_tokenizer.hpp>

TEST(pivot_tokenizer_test, matcher_construction)
{
    using namespace std::literals;
    libio::pivot_matcher m{"pivot"};
    EXPECT_RANGE_EQ(m.needle(), "pivot"s);
}

TEST(pivot_tokenizer_test, matcher_search_full_hit)
{
    using namespace std::literals;
    std::string_view text = "This is a pivotal element!";
    std::span buffer{text.begin(), text.end()};
    libio::pivot_matcher m{"pivot"};
    auto hit = m(buffer);
    EXPECT_RANGE_EQ(hit, "pivot"s);
    EXPECT_EQ(*hit.begin(), 'p');
    EXPECT_EQ(*hit.end(), 'a');
    EXPECT_EQ((hit.begin() - buffer.begin()), 10);
    EXPECT_EQ((buffer.end() - hit.end()), 11);
}

TEST(pivot_tokenizer_test, matcher_search_partial_hit)
{
    using namespace std::literals;
    std::string_view text = "This is an element known as piv";
    std::span buffer{text.begin(), text.end()};
    libio::pivot_matcher m{"pivot"};
    auto hit = m(buffer);
    EXPECT_RANGE_EQ(hit, "piv"s);
    EXPECT_EQ(*hit.begin(), 'p');
    EXPECT_EQ((hit.begin() - buffer.begin()), 28);
    EXPECT_EQ((buffer.end() - hit.end()), 0);
}

TEST(pivot_tokenizer_test, matcher_search_no_hit)
{
    using namespace std::literals;
    std::string_view text = "This has no pivo element wich exists.";
    std::span buffer{text.begin(), text.end()};
    libio::pivot_matcher m{"pivot"};
    auto hit = m(buffer);
    EXPECT_TRUE(std::ranges::empty(hit));
    EXPECT_EQ((hit.begin() - buffer.begin()), text.size());
    EXPECT_EQ((buffer.end() - hit.end()), 0);
}

// Testing the pivot tokenizer
// we need to open a stream with an artifical stream buffer, whose size is set such that we can test overflows.
