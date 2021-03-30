// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <concepts>
#include <type_traits>

 #include <seqan3/test/../../../unit/range/iterator_test_template.hpp>

#include <libjst/utility/bit_vector.hpp>

// ----------------------------------------------------------------------------
// associated types
// ----------------------------------------------------------------------------

TEST(bit_vector_test, value_type)
{
    EXPECT_TRUE((std::same_as<typename libjst::bit_vector<>::value_type, bool>));
}

TEST(bit_vector_test, size_type)
{
    EXPECT_TRUE((std::same_as<typename libjst::bit_vector<>::size_type, size_t>));
}

TEST(bit_vector_test, difference_type)
{
    EXPECT_TRUE((std::integral<typename libjst::bit_vector<>::difference_type>));
}

TEST(bit_vector_test, reference)
{
    EXPECT_TRUE((std::convertible_to<typename libjst::bit_vector<>::reference, bool>));
}

TEST(bit_vector_test, const_reference)
{
    EXPECT_TRUE((std::convertible_to<typename libjst::bit_vector<>::const_reference, bool>));
}

// ----------------------------------------------------------------------------
// Construction and assignment
// ----------------------------------------------------------------------------

TEST(bit_vector_test, construct_with_count)
{
    {
        libjst::bit_vector test_vector(1000);

        EXPECT_EQ(test_vector.size(), 1000u);
    }

    {
        libjst::bit_vector test_vector(64);

        EXPECT_EQ(test_vector.size(), 64u);
    }

    {
        libjst::bit_vector test_vector(512);

        EXPECT_EQ(test_vector.size(), 512u);
    }

    {
        libjst::bit_vector test_vector(1);

        EXPECT_EQ(test_vector.size(), 1u);
    }
}

TEST(bit_vector_test, construct_with_count_and_allocator)
{
    libjst::bit_vector test_vector(1000, true);

    EXPECT_EQ(test_vector.size(), 1000u);

    std::for_each(test_vector.begin(), test_vector.end(), [] (bool const bit)
    {
        EXPECT_TRUE(bit);
    });
}

TEST(bit_vector_test, construct_from_initialiser_list)
{
    { // From the concrete values.
        libjst::bit_vector test_vector{true, false, true, false, false, true, true};
        auto it = test_vector.begin();
        EXPECT_EQ(test_vector.size(), 7u);
        EXPECT_EQ(*it, true);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, true);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, true);
        EXPECT_EQ(*++it, true);
    }

    { // From integers convertible to bool.
        libjst::bit_vector test_vector{1, 0, 1, 0, 0, 1, 1};
        auto it = test_vector.begin();
        EXPECT_EQ(test_vector.size(), 7u);
        EXPECT_EQ(*it, true);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, true);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, false);
        EXPECT_EQ(*++it, true);
        EXPECT_EQ(*++it, true);
    }
}

// ----------------------------------------------------------------------------
// iterators
// ----------------------------------------------------------------------------

TEST(bit_vector_test, begin)
{
    {
        libjst::bit_vector test_vector(1000, true);

        auto it = test_vector.begin();
        EXPECT_TRUE(*it);
    }

    {
        libjst::bit_vector test_vector(64, true);

        auto it = test_vector.begin();
        EXPECT_TRUE(*it);
    }
}

TEST(bit_vector_test, cbegin)
{
    {
        libjst::bit_vector test_vector(1000, true);
        auto cit = std::as_const(test_vector).begin();
        EXPECT_TRUE(*cit);
    }

    {
        libjst::bit_vector test_vector(64, true);
        auto cit = std::ranges::cbegin(std::as_const(test_vector));
        EXPECT_TRUE(*cit);
    }
}

TEST(bit_vector_test, end)
{
    {
        libjst::bit_vector test_vector(1000, true);

        EXPECT_TRUE(test_vector.begin() != test_vector.end());
    }

    {
        libjst::bit_vector test_vector(64, true);
        EXPECT_TRUE(test_vector.begin() != test_vector.end());
    }

    { // empty vector:
        libjst::bit_vector test_vector{};
        EXPECT_TRUE(test_vector.begin() == test_vector.end());
    }
}

TEST(bit_vector_test, cend)
{
    {
        libjst::bit_vector test_vector(1000, true);

        EXPECT_TRUE(std::as_const(test_vector).begin() != std::as_const(test_vector).end());
    }

    {
        libjst::bit_vector test_vector(64, true);
        EXPECT_TRUE(std::as_const(test_vector).begin() != std::ranges::cend(std::as_const(test_vector)));
    }

    { // empty vector:
        libjst::bit_vector test_vector{};
        EXPECT_TRUE(std::as_const(test_vector).begin() == std::as_const(test_vector).end());
    }
}

TEST(bit_vector_test, iterate)
{
    libjst::bit_vector test_vector(70, true);
    auto it = test_vector.begin();
    for (size_t i = 0; i < test_vector.size(); ++i)
    {
        EXPECT_EQ(*it, true);
        ++it;
    }
    EXPECT_TRUE(it == test_vector.end());
}

// ----------------------------------------------------------------------------
// capacity
// ----------------------------------------------------------------------------

TEST(bit_vector_test, size)
{
    {
        libjst::bit_vector test_vector(1000);
        EXPECT_EQ(test_vector.size(), 1000u);
    }

    {
        libjst::bit_vector test_vector(64);
        EXPECT_EQ(test_vector.size(), 64u);
    }

    {
        libjst::bit_vector test_vector(1);
        EXPECT_EQ(test_vector.size(), 1u);
    }

    {
        libjst::bit_vector test_vector{};
        EXPECT_EQ(test_vector.size(), 0u);
    }
}

// ----------------------------------------------------------------------------
// Iterator test
// ----------------------------------------------------------------------------

using bit_vector_iterator = typename libjst::bit_vector<>::iterator;

template <>
struct iterator_fixture<bit_vector_iterator> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;

    static constexpr bool const_iterable = true;

    libjst::bit_vector<> test_range = libjst::bit_vector<>(100, true);
    std::vector<bool> expected_range = std::vector<bool>(100, true);
};

INSTANTIATE_TYPED_TEST_SUITE_P(bit_vector_iterator_test,
                               iterator_fixture,
                               ::testing::Types<bit_vector_iterator>, );

TEST(bit_vector_test, output_iterator)
{
    EXPECT_TRUE((std::output_iterator<typename libjst::bit_vector<>::iterator, bool>));
    EXPECT_FALSE((std::output_iterator<typename libjst::bit_vector<>::const_iterator, bool>));

    libjst::bit_vector test_vector(100, true);
    for (auto it = test_vector.begin(); it != test_vector.end(); it += 2)
        *it = false;

    for (auto it = test_vector.begin(); it != test_vector.end(); ++it)
        EXPECT_EQ(*it, (it - test_vector.begin()) % 2);
}
