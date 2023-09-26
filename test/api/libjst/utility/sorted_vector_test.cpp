// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <type_traits>

#include <libjst/utility/sorted_vector.hpp>

struct sorted_vector_test : public testing::Test
{
    using sorted_vector_t = libjst::sorted_vector<size_t>;
};

// ----------------------------------------------------------------------------
// Modifier
// ----------------------------------------------------------------------------

TEST_F(sorted_vector_test, insert)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.insert(5), 5);
    EXPECT_EQ(*vec.insert(3), 3);
    EXPECT_EQ(*vec.insert(6), 6);
    EXPECT_EQ(*vec.insert(5), 5);
    EXPECT_EQ(*vec.insert(1), 1);
    EXPECT_EQ(*vec.insert(5), 5);
    EXPECT_EQ(*vec.insert(10), 10);
    EXPECT_EQ(*vec.insert(3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
}

TEST_F(sorted_vector_test, insert_hint)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.insert(vec.end(), 5), 5);
    EXPECT_EQ(*vec.insert(vec.begin(), 3), 3);
    EXPECT_EQ(*vec.insert(vec.begin(), 6), 6);
    EXPECT_EQ(*vec.insert(std::ranges::prev(vec.end()), 5), 5);
    EXPECT_EQ(*vec.insert(vec.end(), 1), 1);
    EXPECT_EQ(*vec.insert(vec.begin(), 5), 5);
    EXPECT_EQ(*vec.insert(vec.end(), 10), 10);
    EXPECT_EQ(*vec.insert(std::ranges::next(vec.begin(), 1), 3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
}

TEST_F(sorted_vector_test, emplace)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
}

TEST_F(sorted_vector_test, emplace_hint)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace_hint(vec.end(), 5), 5);
    EXPECT_EQ(*vec.emplace_hint(vec.begin(), 3), 3);
    EXPECT_EQ(*vec.emplace_hint(vec.begin(), 6), 6);
    EXPECT_EQ(*vec.emplace_hint(std::ranges::prev(vec.end()), 5), 5);
    EXPECT_EQ(*vec.emplace_hint(vec.end(), 1), 1);
    EXPECT_EQ(*vec.emplace_hint(vec.begin(), 5), 5);
    EXPECT_EQ(*vec.emplace_hint(vec.end(), 10), 10);
    EXPECT_EQ(*vec.emplace_hint(std::ranges::next(vec.begin(), 1), 3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
}

TEST_F(sorted_vector_test, erase)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
    EXPECT_EQ(*vec.erase(vec.begin()), 3);
    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{3, 3, 5, 5, 5, 6, 10})));
    EXPECT_EQ(*vec.erase(std::ranges::next(vec.begin())), 5);
    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{3, 5, 5, 5, 6, 10})));
    auto it = vec.erase(std::ranges::next(vec.begin(), 5));
    EXPECT_TRUE(it == vec.end());
    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{3, 5, 5, 5, 6})));
}

TEST_F(sorted_vector_test, erase_range)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 3, 3, 5, 5, 5, 6, 10})));
    EXPECT_EQ(*vec.erase(std::ranges::next(vec.begin()), std::ranges::next(vec.begin(), 6)), 6);
    EXPECT_TRUE(std::ranges::equal(vec, (std::vector<size_t>{1, 6, 10})));

    vec.erase(vec.begin(), vec.end());
    EXPECT_TRUE(vec.empty());
}

TEST_F(sorted_vector_test, clear)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_FALSE(vec.empty());
    vec.clear();
    EXPECT_TRUE(vec.empty());
}

// ----------------------------------------------------------------------------
// Cpacity
// ----------------------------------------------------------------------------

TEST_F(sorted_vector_test, empty)
{
    sorted_vector_t vec{};
    EXPECT_TRUE(vec.empty());

    vec.insert(5);
    EXPECT_FALSE(std::as_const(vec).empty());
}

TEST_F(sorted_vector_test, size)
{
    sorted_vector_t vec{};
    EXPECT_EQ(vec.size(), 0u);

    vec.insert(5);
    EXPECT_EQ(std::as_const(vec).size(), 1u);
    vec.insert(5);
    vec.insert(6);
    vec.insert(1);
    EXPECT_EQ(vec.size(), 4u);
}

TEST_F(sorted_vector_test, max_size)
{
    sorted_vector_t vec{};
    EXPECT_EQ(std::as_const(vec).max_size(), std::vector<size_t>{}.max_size());
}

// ----------------------------------------------------------------------------
// Lookup
// ----------------------------------------------------------------------------

TEST_F(sorted_vector_test, find)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_EQ(*vec.find(5), 5);
    EXPECT_EQ(*std::as_const(vec).find(6), 6);
    EXPECT_TRUE(vec.find(7) == vec.end());
}

TEST_F(sorted_vector_test, contains)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_TRUE(vec.contains(5));
    EXPECT_TRUE(std::as_const(vec).contains(1));
    EXPECT_TRUE(std::as_const(vec).contains(6));
    EXPECT_TRUE(std::as_const(vec).contains(3));
    EXPECT_FALSE(std::as_const(vec).contains(7));
    EXPECT_FALSE(std::as_const(vec).contains(11));
    EXPECT_FALSE(std::as_const(vec).contains(0));
}

TEST_F(sorted_vector_test, equal_range)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    { // no key
        auto [first, last] = vec.equal_range(0);
        EXPECT_EQ(std::ranges::distance(first, last), 0);
    }

    { // key = 1
        auto [first, last] = std::as_const(vec).equal_range(1);
        EXPECT_EQ(*first, 1);
        EXPECT_EQ(std::ranges::distance(first, last), 1);
    }

    { // key = 3
        auto [first, last] = std::as_const(vec).equal_range(3);
        EXPECT_EQ(*first, 3);
        EXPECT_EQ(*std::ranges::next(first), 3);
        EXPECT_EQ(std::ranges::distance(first, last), 2);
    }

    { // key = 5
        auto [first, last] = std::as_const(vec).equal_range(5);
        EXPECT_EQ(*first, 5);
        EXPECT_EQ(*std::ranges::next(first), 5);
        EXPECT_EQ(*std::ranges::next(first, 2), 5);
        EXPECT_EQ(std::ranges::distance(first, last), 3);
    }

    { // key = 6
        auto [first, last] = std::as_const(vec).equal_range(6);
        EXPECT_EQ(*first, 6);
        EXPECT_EQ(std::ranges::distance(first, last), 1);
    }

    { // key = 7
        auto [first, last] = std::as_const(vec).equal_range(7);
        EXPECT_EQ(std::ranges::distance(first, last), 0);
    }

    { // key = 10
        auto [first, last] = std::as_const(vec).equal_range(10);
        EXPECT_EQ(*first, 10);
        EXPECT_EQ(std::ranges::distance(first, last), 1);
    }

    { // key = 11
        auto [first, last] = std::as_const(vec).equal_range(11);
        EXPECT_EQ(std::ranges::distance(first, last), 0);
    }
}

TEST_F(sorted_vector_test, count)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_EQ(vec.count(0), 0u);
    EXPECT_EQ(std::as_const(vec).count(1), 1u);
    EXPECT_EQ(std::as_const(vec).count(2), 0u);
    EXPECT_EQ(std::as_const(vec).count(3), 2u);
    EXPECT_EQ(std::as_const(vec).count(5), 3u);
    EXPECT_EQ(std::as_const(vec).count(10), 1u);
    EXPECT_EQ(std::as_const(vec).count(11), 0u);
}

TEST_F(sorted_vector_test, lower_bound)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    EXPECT_EQ(*vec.lower_bound(0), 1u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(1), 1u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(2), 3u);
    auto it = std::as_const(vec).lower_bound(3);
    EXPECT_EQ(*it, 3u);
    EXPECT_EQ(*std::ranges::next(it), 3u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(4), 5u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(5), 5u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(6), 6u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(7), 10u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(8), 10u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(9), 10u);
    EXPECT_EQ(*std::as_const(vec).lower_bound(10), 10u);
    EXPECT_TRUE(std::as_const(vec).lower_bound(11) == vec.end());
}

TEST_F(sorted_vector_test, upper_bound)
{
    sorted_vector_t vec{};

    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(3), 3);
    EXPECT_EQ(*vec.emplace(6), 6);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(1), 1);
    EXPECT_EQ(*vec.emplace(5), 5);
    EXPECT_EQ(*vec.emplace(10), 10);
    EXPECT_EQ(*vec.emplace(3), 3);

    // first element that is greater than x
    EXPECT_EQ(*vec.upper_bound(0), 1u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(1), 3u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(2), 3u);
    auto it = std::as_const(vec).upper_bound(3);
    EXPECT_EQ(*it, 5u);
    EXPECT_EQ(*++it, 5u);
    EXPECT_EQ(*++it, 5u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(4), 5u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(5), 6u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(6), 10u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(7), 10u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(8), 10u);
    EXPECT_EQ(*std::as_const(vec).upper_bound(9), 10u);
    EXPECT_TRUE(std::as_const(vec).upper_bound(10) == vec.end());
    EXPECT_TRUE(std::as_const(vec).upper_bound(11) == vec.end());
}
