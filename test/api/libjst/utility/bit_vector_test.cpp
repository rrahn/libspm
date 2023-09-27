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

template <typename t>
struct bit_vector_test : public testing::Test
{
    using bit_vector_type = t;
};

// Next, associate a list of types with the test suite, which will be
// repeated for each type in the list.  The typedef is necessary for
// the macro to parse correctly.
using bit_vector_types = testing::Types<libjst::bit_vector<>>;
TYPED_TEST_SUITE(bit_vector_test, bit_vector_types, );

// ----------------------------------------------------------------------------
// associated types
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, value_type)
{
    EXPECT_TRUE((std::same_as<typename TypeParam::value_type, bool>));
}

TYPED_TEST(bit_vector_test, size_type)
{
    EXPECT_TRUE((std::same_as<typename TypeParam::size_type, size_t>));
}

TYPED_TEST(bit_vector_test, difference_type)
{
    EXPECT_TRUE((std::integral<typename TypeParam::difference_type>));
}

TYPED_TEST(bit_vector_test, reference)
{
    EXPECT_TRUE((std::convertible_to<typename TypeParam::reference, bool>));
}

TYPED_TEST(bit_vector_test, const_reference)
{
    EXPECT_TRUE((std::convertible_to<typename TypeParam::const_reference, bool>));
}

// ----------------------------------------------------------------------------
// Construction and assignment
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, construct_with_count)
{
    {
        TypeParam test_vector(1000);

        EXPECT_EQ(test_vector.size(), 1000u);
    }

    {
        TypeParam test_vector(64);

        EXPECT_EQ(test_vector.size(), 64u);
    }

    {
        TypeParam test_vector(512);

        EXPECT_EQ(test_vector.size(), 512u);
    }

    {
        TypeParam test_vector(1);

        EXPECT_EQ(test_vector.size(), 1u);
    }
}

TYPED_TEST(bit_vector_test, construct_with_count_and_allocator)
{
    TypeParam test_vector(1000, true);

    EXPECT_EQ(test_vector.size(), 1000u);

    std::for_each(test_vector.begin(), test_vector.end(), [] (bool const bit)
    {
        EXPECT_TRUE(bit);
    });
}

TYPED_TEST(bit_vector_test, construct_from_initialiser_list)
{
    { // From the concrete values.
        TypeParam test_vector{true, false, true, false, false, true, true};
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
        TypeParam test_vector{1, 0, 1, 0, 0, 1, 1};
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

TYPED_TEST(bit_vector_test, assign_with_count)
{
    {
        TypeParam test_vector{};
        test_vector.assign(1000, false);

        EXPECT_EQ(test_vector.size(), 1000u);
    }

    {
        TypeParam test_vector{};
        test_vector.assign(64, true);

        EXPECT_EQ(test_vector.size(), 64u);
    }

    {
        TypeParam test_vector{};
        test_vector.assign(0, false);

        EXPECT_EQ(test_vector.size(), 0u);
    }

    {
        TypeParam test_vector{};
        test_vector.assign(1, true);

        EXPECT_EQ(test_vector.size(), 1u);
    }
}

TYPED_TEST(bit_vector_test, assign_from_initialiser_list)
{
    { // From the concrete values.
        TypeParam test_vector{};
        test_vector.assign({true, false, true, false, false, true, true});

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
        TypeParam test_vector{};
        test_vector.assign({1, 0, 1, 0, 0, 1, 1});

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

    {
        TypeParam test_vector{};
        test_vector.assign(std::initializer_list<bool>{});

        EXPECT_EQ(test_vector.size(), 0u);
    }
}

// ----------------------------------------------------------------------------
// iterators
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, begin)
{
    {
        TypeParam test_vector(1000, true);

        auto it = test_vector.begin();
        EXPECT_TRUE(*it);
    }

    {
        TypeParam test_vector(64, true);

        auto it = test_vector.begin();
        EXPECT_TRUE(*it);
    }
}

TYPED_TEST(bit_vector_test, cbegin)
{
    {
        TypeParam test_vector(1000, true);
        auto cit = std::as_const(test_vector).begin();
        EXPECT_TRUE(*cit);
    }

    {
        TypeParam test_vector(64, true);
        auto cit = std::ranges::cbegin(std::as_const(test_vector));
        EXPECT_TRUE(*cit);
    }
}

TYPED_TEST(bit_vector_test, end)
{
    {
        TypeParam test_vector(1000, true);

        EXPECT_TRUE(test_vector.begin() != test_vector.end());
    }

    {
        TypeParam test_vector(64, true);
        EXPECT_TRUE(test_vector.begin() != test_vector.end());
    }

    { // empty vector:
        TypeParam test_vector{};
        EXPECT_TRUE(test_vector.begin() == test_vector.end());
    }
}

TYPED_TEST(bit_vector_test, cend)
{
    {
        TypeParam test_vector(1000, true);

        EXPECT_TRUE(std::as_const(test_vector).begin() != std::as_const(test_vector).end());
    }

    {
        TypeParam test_vector(64, true);
        EXPECT_TRUE(std::as_const(test_vector).begin() != std::ranges::cend(std::as_const(test_vector)));
    }

    { // empty vector:
        TypeParam test_vector{};
        EXPECT_TRUE(std::as_const(test_vector).begin() == std::as_const(test_vector).end());
    }
}

TYPED_TEST(bit_vector_test, iterate)
{
    TypeParam test_vector(70, true);
    auto it = test_vector.begin();
    for (size_t i = 0; i < test_vector.size(); ++i)
    {
        EXPECT_EQ(*it, true);
        ++it;
    }
    EXPECT_TRUE(it == test_vector.end());
}

// ----------------------------------------------------------------------------
// Element access
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, subscript_operator)
{
    TypeParam test_vector{true, false, true, false, false, true, true};
    EXPECT_EQ(test_vector.size(), 7u);
    EXPECT_EQ(test_vector[0], true);
    EXPECT_EQ(test_vector[1], false);
    EXPECT_EQ(test_vector[2], true);
    EXPECT_EQ(std::as_const(test_vector)[3], false);
    EXPECT_EQ(std::as_const(test_vector)[4], false);
    EXPECT_EQ(std::as_const(test_vector)[5], true);
    EXPECT_EQ(std::as_const(test_vector)[6], true);
}

TYPED_TEST(bit_vector_test, back)
{
    TypeParam test_vector{true, false, true, false, false, true, true};
    EXPECT_EQ(test_vector.size(), 7u);
    EXPECT_EQ(test_vector.back(), true);
    EXPECT_EQ(std::as_const(test_vector).back(), true);
}

TYPED_TEST(bit_vector_test, all)
{
    { // empty vector
        TypeParam test_vector{};
        EXPECT_TRUE(test_vector.all());
    }

    {
        TypeParam test_vector(250, true);
        EXPECT_TRUE(test_vector.all());

        test_vector[249] = false;
        EXPECT_FALSE(test_vector.all());

        test_vector[249] = true;
        test_vector[0] = false;
        EXPECT_FALSE(test_vector.all());
    }
}

TYPED_TEST(bit_vector_test, any)
{
    { // empty vector
        TypeParam test_vector{};
        EXPECT_FALSE(test_vector.any());
    }

    {
        TypeParam test_vector(250, false);
        EXPECT_FALSE(test_vector.any());

        test_vector[249] = true;
        EXPECT_TRUE(test_vector.any());

        test_vector[249] = false;
        test_vector[0] = true;
        EXPECT_TRUE(test_vector.any());
    }
}

TYPED_TEST(bit_vector_test, none)
{
    { // empty vector
        TypeParam test_vector{};
        EXPECT_TRUE(test_vector.none());
    }

    {
        TypeParam test_vector(250, false);
        EXPECT_TRUE(test_vector.none());

        test_vector[249] = true;
        EXPECT_FALSE(test_vector.none());

        test_vector[249] = false;
        test_vector[0] = true;
        EXPECT_FALSE(test_vector.none());
    }
}

// ----------------------------------------------------------------------------
// Modifiers
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, resize)
{
    TypeParam test_vector{};

    EXPECT_EQ(test_vector.size(), 0u);
    test_vector.resize(64);
    EXPECT_EQ(test_vector.size(), 64u);
    EXPECT_TRUE(test_vector.none());

    test_vector.resize(128, true);
    EXPECT_EQ(test_vector.size(), 128u);
    EXPECT_TRUE(test_vector.any());

    test_vector.resize(1, true);
    EXPECT_EQ(test_vector.size(), 1u);
    EXPECT_TRUE(test_vector.none());

    test_vector.resize(0, true);
    EXPECT_EQ(test_vector.size(), 0u);
    EXPECT_TRUE(test_vector.none());
}

TYPED_TEST(bit_vector_test, push_back)
{
    TypeParam test_vector{};

    EXPECT_EQ(test_vector.size(), 0u);
    test_vector.push_back(true);
    EXPECT_EQ(test_vector.size(), 1u);
    EXPECT_TRUE(test_vector.back());

    test_vector.resize(128, true);
    test_vector.push_back(false);
    EXPECT_EQ(test_vector.size(), 129u);
    EXPECT_FALSE(test_vector.back());
}

TYPED_TEST(bit_vector_test, swap)
{
    TypeParam test_vector_left{};
    TypeParam test_vector_right(250, true);

    test_vector_left.swap(test_vector_right);

    EXPECT_EQ(test_vector_left.size(), 250u);
    EXPECT_EQ(test_vector_right.size(), 0u);
    EXPECT_TRUE(test_vector_left.all());

    test_vector_right.resize(78);
    test_vector_left.swap(test_vector_right);
    EXPECT_EQ(test_vector_left.size(), 78u);
    EXPECT_EQ(test_vector_right.size(), 250u);
    EXPECT_TRUE(test_vector_left.none());
    EXPECT_TRUE(test_vector_right.all());
}

TYPED_TEST(bit_vector_test, operator_binary_and)
{
    TypeParam test_vector(250, false);
    TypeParam test_vector_all(250, true);

    test_vector &= test_vector_all;

    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_EQ(test_vector_all.size(), 250u);

    EXPECT_TRUE(test_vector.none());
    EXPECT_TRUE(test_vector_all.all());

    test_vector[0] = true;
    test_vector[10] = true;
    test_vector[63] = true;
    test_vector[64] = true;
    test_vector[127] = true;
    test_vector[128] = true;
    test_vector[200] = true;
    test_vector[249] = true;

    test_vector = test_vector & test_vector_all;
    EXPECT_FALSE(test_vector.none());
    EXPECT_TRUE(test_vector_all.all());

    EXPECT_TRUE(test_vector[0]);
    EXPECT_FALSE(test_vector[1]);
    EXPECT_FALSE(test_vector[9]);
    EXPECT_TRUE(test_vector[10]);
    EXPECT_FALSE(test_vector[11]);
    EXPECT_FALSE(test_vector[62]);
    EXPECT_TRUE(test_vector[63]);
    EXPECT_TRUE(test_vector[64]);
    EXPECT_FALSE(test_vector[65]);
    EXPECT_FALSE(test_vector[126]);
    EXPECT_TRUE(test_vector[127]);
    EXPECT_TRUE(test_vector[128]);
    EXPECT_FALSE(test_vector[129]);
    EXPECT_FALSE(test_vector[199]);
    EXPECT_TRUE(test_vector[200]);
    EXPECT_FALSE(test_vector[201]);
    EXPECT_FALSE(test_vector[248]);
    EXPECT_TRUE(test_vector[249]);
}

TYPED_TEST(bit_vector_test, operator_binary_or)
{
    TypeParam test_vector(250, false);
    TypeParam test_vector_all(250, true);

    test_vector |= test_vector_all;

    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_EQ(test_vector_all.size(), 250u);

    EXPECT_TRUE(test_vector.all());
    EXPECT_TRUE(test_vector_all.all());

    test_vector[0] = false;
    test_vector[10] = false;
    test_vector[63] = false;
    test_vector[64] = false;
    test_vector[127] = false;
    test_vector[128] = false;
    test_vector[200] = false;
    test_vector[249] = false;

    test_vector_all[0] = false;
    test_vector_all[10] = false;
    test_vector_all[127] = false;
    test_vector_all[128] = false;
    test_vector_all[249] = false;

    test_vector = test_vector | test_vector_all;
    EXPECT_FALSE(test_vector.all());
    EXPECT_FALSE(test_vector_all.all());

    EXPECT_FALSE(test_vector[0]);
    EXPECT_TRUE(test_vector[1]);
    EXPECT_TRUE(test_vector[9]);
    EXPECT_FALSE(test_vector[10]);
    EXPECT_TRUE(test_vector[11]);
    EXPECT_TRUE(test_vector[62]);
    EXPECT_TRUE(test_vector[63]);
    EXPECT_TRUE(test_vector[64]);
    EXPECT_TRUE(test_vector[65]);
    EXPECT_TRUE(test_vector[126]);
    EXPECT_FALSE(test_vector[127]);
    EXPECT_FALSE(test_vector[128]);
    EXPECT_TRUE(test_vector[129]);
    EXPECT_TRUE(test_vector[199]);
    EXPECT_TRUE(test_vector[200]);
    EXPECT_TRUE(test_vector[201]);
    EXPECT_TRUE(test_vector[248]);
    EXPECT_FALSE(test_vector[249]);
}

TYPED_TEST(bit_vector_test, operator_binary_xor)
{
    TypeParam test_vector(250, false);
    TypeParam test_vector_all(250, true);

    test_vector ^= test_vector_all;

    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_EQ(test_vector_all.size(), 250u);

    EXPECT_TRUE(test_vector.all());
    EXPECT_TRUE(test_vector_all.all());

    test_vector[0] = false;
    test_vector[10] = false;
    test_vector[63] = false;
    test_vector[64] = false;
    test_vector[127] = false;
    test_vector[128] = false;
    test_vector[200] = false;
    test_vector[249] = false;

    test_vector_all[0] = false;
    test_vector_all[10] = false;
    test_vector_all[127] = false;
    test_vector_all[128] = false;
    test_vector_all[249] = false;

    test_vector = test_vector ^ test_vector_all;
    EXPECT_FALSE(test_vector.all());
    EXPECT_FALSE(test_vector_all.all());

    EXPECT_FALSE(test_vector[0]);
    EXPECT_FALSE(test_vector[1]);
    EXPECT_FALSE(test_vector[9]);
    EXPECT_FALSE(test_vector[10]);
    EXPECT_FALSE(test_vector[11]);
    EXPECT_FALSE(test_vector[62]);
    EXPECT_TRUE(test_vector[63]);
    EXPECT_TRUE(test_vector[64]);
    EXPECT_FALSE(test_vector[65]);
    EXPECT_FALSE(test_vector[126]);
    EXPECT_FALSE(test_vector[127]);
    EXPECT_FALSE(test_vector[128]);
    EXPECT_FALSE(test_vector[129]);
    EXPECT_FALSE(test_vector[199]);
    EXPECT_TRUE(test_vector[200]);
    EXPECT_FALSE(test_vector[201]);
    EXPECT_FALSE(test_vector[248]);
    EXPECT_FALSE(test_vector[249]);
}

TYPED_TEST(bit_vector_test, operator_binary_not)
{
    TypeParam test_vector(250, false);
    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_TRUE(test_vector.none());

    TypeParam expected_vector = ~test_vector;
    EXPECT_EQ(expected_vector.size(), 250u);
    EXPECT_TRUE(test_vector.none());
    EXPECT_TRUE(expected_vector.all());

    test_vector[0] = true;
    test_vector[10] = true;
    test_vector[63] = true;
    test_vector[64] = true;
    test_vector[127] = true;
    test_vector[128] = true;
    test_vector[200] = true;
    test_vector[249] = true;

    expected_vector = ~test_vector;
    EXPECT_FALSE(expected_vector.all());

    EXPECT_FALSE(expected_vector[0]);
    EXPECT_TRUE(expected_vector[1]);
    EXPECT_TRUE(expected_vector[9]);
    EXPECT_FALSE(expected_vector[10]);
    EXPECT_TRUE(expected_vector[11]);
    EXPECT_TRUE(expected_vector[62]);
    EXPECT_FALSE(expected_vector[63]);
    EXPECT_FALSE(expected_vector[64]);
    EXPECT_TRUE(expected_vector[65]);
    EXPECT_TRUE(expected_vector[126]);
    EXPECT_FALSE(expected_vector[127]);
    EXPECT_FALSE(expected_vector[128]);
    EXPECT_TRUE(expected_vector[129]);
    EXPECT_TRUE(expected_vector[199]);
    EXPECT_FALSE(expected_vector[200]);
    EXPECT_TRUE(expected_vector[201]);
    EXPECT_TRUE(expected_vector[248]);
    EXPECT_FALSE(expected_vector[249]);
}

TYPED_TEST(bit_vector_test, flip)
{
    TypeParam test_vector(250, false);
    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_TRUE(test_vector.none());

    test_vector.flip();
    EXPECT_TRUE(test_vector.all());

    test_vector[0] = false;
    test_vector[10] = false;
    test_vector[63] = false;
    test_vector[64] = false;
    test_vector[127] = false;
    test_vector[128] = false;
    test_vector[200] = false;
    test_vector[249] = false;

    test_vector.flip();

    EXPECT_TRUE(test_vector[0]);
    EXPECT_FALSE(test_vector[1]);
    EXPECT_FALSE(test_vector[9]);
    EXPECT_TRUE(test_vector[10]);
    EXPECT_FALSE(test_vector[11]);
    EXPECT_FALSE(test_vector[62]);
    EXPECT_TRUE(test_vector[63]);
    EXPECT_TRUE(test_vector[64]);
    EXPECT_FALSE(test_vector[65]);
    EXPECT_FALSE(test_vector[126]);
    EXPECT_TRUE(test_vector[127]);
    EXPECT_TRUE(test_vector[128]);
    EXPECT_FALSE(test_vector[129]);
    EXPECT_FALSE(test_vector[199]);
    EXPECT_TRUE(test_vector[200]);
    EXPECT_FALSE(test_vector[201]);
    EXPECT_FALSE(test_vector[248]);
    EXPECT_TRUE(test_vector[249]);
}

TYPED_TEST(bit_vector_test, flip_single_bit)
{
    TypeParam test_vector(250, false);
    EXPECT_EQ(test_vector.size(), 250u);
    EXPECT_TRUE(test_vector.none());

    test_vector.flip(0);
    test_vector.flip(10);
    test_vector.flip(63);
    test_vector.flip(64);
    test_vector.flip(127);
    test_vector.flip(128);
    test_vector.flip(200);
    test_vector.flip(249);

    EXPECT_TRUE(test_vector[0]);
    EXPECT_FALSE(test_vector[1]);
    EXPECT_FALSE(test_vector[9]);
    EXPECT_TRUE(test_vector[10]);
    EXPECT_FALSE(test_vector[11]);
    EXPECT_FALSE(test_vector[62]);
    EXPECT_TRUE(test_vector[63]);
    EXPECT_TRUE(test_vector[64]);
    EXPECT_FALSE(test_vector[65]);
    EXPECT_FALSE(test_vector[126]);
    EXPECT_TRUE(test_vector[127]);
    EXPECT_TRUE(test_vector[128]);
    EXPECT_FALSE(test_vector[129]);
    EXPECT_FALSE(test_vector[199]);
    EXPECT_TRUE(test_vector[200]);
    EXPECT_FALSE(test_vector[201]);
    EXPECT_FALSE(test_vector[248]);
    EXPECT_TRUE(test_vector[249]);

    test_vector.flip(0);
    test_vector.flip(10);
    test_vector.flip(63);
    test_vector.flip(64);

    EXPECT_FALSE(test_vector[0]);
    EXPECT_FALSE(test_vector[10]);
    EXPECT_FALSE(test_vector[63]);
    EXPECT_FALSE(test_vector[64]);

    EXPECT_THROW(test_vector.flip(250), std::out_of_range);
}

// ----------------------------------------------------------------------------
// capacity
// ----------------------------------------------------------------------------

TYPED_TEST(bit_vector_test, size)
{
    {
        TypeParam test_vector(1000);
        EXPECT_EQ(test_vector.size(), 1000u);
    }

    {
        TypeParam test_vector(64);
        EXPECT_EQ(test_vector.size(), 64u);
    }

    {
        TypeParam test_vector(1);
        EXPECT_EQ(test_vector.size(), 1u);
    }

    {
        TypeParam test_vector{};
        EXPECT_EQ(test_vector.size(), 0u);
    }
}

TYPED_TEST(bit_vector_test, empty)
{
    {
        TypeParam test_vector{};
        EXPECT_TRUE(test_vector.empty());
    }

    {
        TypeParam test_vector(1);
        EXPECT_FALSE(test_vector.empty());
    }

    {
        TypeParam test_vector(1000);
        EXPECT_FALSE(test_vector.empty());
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

TYPED_TEST(bit_vector_test, output_iterator)
{
    EXPECT_TRUE((std::output_iterator<typename libjst::bit_vector<>::iterator, bool>));
    EXPECT_FALSE((std::output_iterator<typename libjst::bit_vector<>::const_iterator, bool>));

    libjst::bit_vector<> test_vector(100, true);
    for (auto it = test_vector.begin(); it != test_vector.end(); it += 2)
        *it = false;

    for (auto it = test_vector.begin(); it != test_vector.end(); ++it)
        EXPECT_EQ(*it, (it - test_vector.begin()) % 2);
}
