// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>
#include <string>
#include <type_traits>

#include <seqan3/alphabet/adaptation/char.hpp>

#include <libjst/container/journaled_sequence_group.hpp>
#include <libjst/variant/variant_generic.hpp>
#include <libjst/variant/variant_store_covered.hpp>
#include <libjst/utility/bit_vector.hpp>

template <typename t>
struct journaled_sequence_group_fixture : public ::testing::Test
{
    using source_t = typename t::first_type;
    using store_t = typename t::second_type;

    using type = libjst::journaled_sequence_group<std::views::all_t<source_t const &>, store_t>;

    source_t _test_source{"test source sequence"};

    source_t const & source() const noexcept
    {
        return _test_source;
    }
};

using test_types = ::testing::Types<
    std::pair<std::string,
              libjst::variant_store_covered<std::vector<libjst::generic_variant<char>>, libjst::bit_vector<>>>
>;
TYPED_TEST_SUITE(journaled_sequence_group_fixture, test_types, );

TYPED_TEST(journaled_sequence_group_fixture, constructibility)
{
    using type = typename TestFixture::type;

    ASSERT_TRUE(std::is_default_constructible_v<type>);
    ASSERT_TRUE(std::is_copy_constructible_v<type>);
    ASSERT_TRUE(std::is_nothrow_move_constructible_v<type>);
    ASSERT_TRUE(std::is_copy_assignable_v<type>);
    ASSERT_TRUE(std::is_nothrow_move_assignable_v<type>);
    ASSERT_TRUE(std::is_destructible_v<type>);
}

TYPED_TEST(journaled_sequence_group_fixture, construction)
{
    FAIL();
    // From source sequence and number of journaled sequences (empty by default?)
    // ASSERT_TRUE((std::is_constructible_v<type, std::views::all_t<source_t>, size_t>));
    // // From source sequence and predefined set of variants (is owned later by group).
    // ASSERT_TRUE((std::is_constructible_v<type, std::views::all_t<source_t>, store_t const &>));
    // ASSERT_FALSE((std::is_constructible_v<type, std::string, store_t const &>)); // source must be viewable_range
}

TYPED_TEST(journaled_sequence_group_fixture, clear)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, reset)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, at)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, size)
{
    using type = typename TestFixture::type;
    type js_group{};

    EXPECT_EQ(js_group.size(), 0u);
    js_group = type{this->source(), 10u};
    EXPECT_EQ(js_group.size(), 10u);
}

TYPED_TEST(journaled_sequence_group_fixture, load)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, save)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, begin)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, end)
{
    FAIL();
}
