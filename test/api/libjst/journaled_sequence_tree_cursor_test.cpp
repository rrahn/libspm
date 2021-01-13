// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp> // allow std::string be recognised as seqan3::sequence
#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/journaled_sequence_tree.hpp>
#include <libjst/journaled_sequence_tree_cursor.hpp>

using namespace std::literals;

struct journaled_sequence_tree_cursor_fixture : public ::testing::Test
{
    using sequence_t = std::string;
    using jst_t = libjst::journaled_sequence_tree<sequence_t>;
    using jst_cursor_t = libjst::journaled_sequence_tree_cursor<jst_t>;
    using context_position_t = typename jst_cursor_t::context_position_type;

    using aligned_sequence_t = std::vector<seqan3::gapped<char>>;
    using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;

    sequence_t reference{"aaaabbbbcccc"};

    alignment_t alignment1{make_gapped("aaaabbbbcccc------"sv), make_gapped("------------aabbcc"sv)};
    alignment_t alignment2{make_gapped("aaaabbbbcccc------"sv), make_gapped("------------abcabc"sv)};
    alignment_t alignment3{make_gapped("aaaa--bbbb--cccc--"sv), make_gapped("----cc----aa----bb"sv)};

    jst_t jst{};

    void SetUp() override
    {
        jst = jst_t{sequence_t{reference}};
        jst.add(alignment1);
        jst.add(alignment2);
        jst.add(alignment3);
    }

    static aligned_sequence_t make_gapped(std::string_view const seq)
    {
        aligned_sequence_t tmp{};
        tmp.reserve(seq.size());

        std::for_each(seq.begin(), seq.end(), [&] (char const c)
        {
            if (c == '-')
                tmp.emplace_back(seqan3::gap{});
            else
                tmp.emplace_back(c);
        });

        return tmp;
    }
};

TEST_F(journaled_sequence_tree_cursor_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<jst_cursor_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_cursor_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<jst_cursor_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_cursor_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<jst_cursor_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<jst_cursor_t>);
    EXPECT_TRUE((std::is_nothrow_constructible_v<jst_cursor_t, jst_t const *, size_t const>));
}

TEST_F(journaled_sequence_tree_cursor_fixture, context)
{
    jst_cursor_t jst_cursor{std::addressof(jst), 4u};

    EXPECT_RANGE_EQ(jst_cursor.context(), "aabb"sv);
    EXPECT_RANGE_EQ(std::as_const(jst_cursor).context(), "aabb"sv);
}

TEST_F(journaled_sequence_tree_cursor_fixture, positions)
{
    jst_cursor_t jst_cursor{std::addressof(jst), 4u};
    EXPECT_RANGE_EQ(jst_cursor.positions(), (std::vector{context_position_t{0u, 0u}}));
    EXPECT_RANGE_EQ(std::as_const(jst_cursor).positions(), (std::vector{context_position_t{0u, 0u}}));
}

TEST_F(journaled_sequence_tree_cursor_fixture, advance)
{
    jst_cursor_t jst_cursor{std::addressof(jst), 4u};

    EXPECT_RANGE_EQ(std::as_const(jst_cursor).context(), "aabb"sv);
    EXPECT_RANGE_EQ(std::as_const(jst_cursor).positions(), (std::vector{context_position_t{0u, 0u}}));

    jst_cursor.advance();

    EXPECT_RANGE_EQ(std::as_const(jst_cursor).context(), "abbc"sv);
    EXPECT_RANGE_EQ(std::as_const(jst_cursor).positions(), (std::vector{context_position_t{0u, 1u}}));
}

TEST_F(journaled_sequence_tree_cursor_fixture, at_end)
{
    jst_cursor_t jst_cursor{std::addressof(jst), 4u};

    EXPECT_FALSE(jst_cursor.at_end());

    for (unsigned i = 0; i < 8; ++i)
    {
        jst_cursor.advance();
        EXPECT_FALSE(jst_cursor.at_end());
    }

    jst_cursor.advance();
    EXPECT_TRUE(std::as_const(jst_cursor).at_end());
}

TEST_F(journaled_sequence_tree_cursor_fixture, context_empty_jst)
{
    jst_t jst_empty{};
    jst_cursor_t jst_cursor{std::addressof(jst_empty), 4u};

    EXPECT_TRUE(jst_cursor.at_end());
}

TEST_F(journaled_sequence_tree_cursor_fixture, context_too_large)
{
    jst_cursor_t jst_cursor{std::addressof(jst), 7u};

    EXPECT_TRUE(jst_cursor.at_end());
}
