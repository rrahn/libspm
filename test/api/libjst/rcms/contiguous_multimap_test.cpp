// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/rcms/contiguous_multimap.hpp>

struct contiguous_multimap_test : public ::testing::Test {
    using key_type = uint32_t;
    using mapped_type = uint32_t;
    using map_type = libjst::contiguous_multimap<key_type, mapped_type>;
    using value_type = std::ranges::range_value_t<map_type>;
};

TEST_F(contiguous_multimap_test, insert) {
    map_type map{};
    auto it = map.insert(value_type{10, 0});
    EXPECT_EQ(it->second, 0);
    EXPECT_EQ(it->first, 10);

    it = map.insert(value_type{25, 1});
    EXPECT_EQ(it->second, 1);
    EXPECT_EQ(it->first, 25);

    it = map.insert(value_type{3, 3});
    EXPECT_EQ(it->second, 3);
    EXPECT_EQ(it->first, 3);

    it = map.insert(value_type{25, 2});
    EXPECT_EQ(it->second, 2);
    EXPECT_EQ(it->first, 25);
}

TEST_F(contiguous_multimap_test, iterate) {
    map_type map{};

    map.insert(value_type{10, 0});
    map.insert(value_type{25, 1});
    map.insert(value_type{3, 3});
    map.insert(value_type{25, 2});

    auto it = map.begin();
    EXPECT_EQ(it->first, 3);
    EXPECT_EQ(it->second, 3);
    ++it;
    EXPECT_EQ(it->first, 10);
    EXPECT_EQ(it->second, 0);
    ++it;
    EXPECT_EQ(it->first, 25);
    EXPECT_EQ(it->second, 1);
    ++it;
    EXPECT_EQ(it->first, 25);
    EXPECT_EQ(it->second, 2);
    ++it;
    EXPECT_TRUE(it == map.end());
}

TEST_F(contiguous_multimap_test, reference) {
    map_type map{};

    map.insert(value_type{10, 0});
    map.insert(value_type{25, 1});
    map.insert(value_type{3, 3});
    map.insert(value_type{25, 2});

    auto it = map.begin();
    it->second = 1;
    EXPECT_EQ(it->first, 3);
    EXPECT_EQ(it->second, 1);
    ++it;
    it->second = 3;
    EXPECT_EQ(it->first, 10);
    EXPECT_EQ(it->second, 3);
    ++it;
    it->second = 0;
    EXPECT_EQ(it->first, 25);
    EXPECT_EQ(it->second, 0);
    ++it;
    it->second = 2;
    EXPECT_EQ(it->first, 25);
    EXPECT_EQ(it->second, 2);
    EXPECT_TRUE(++it == map.end());
}

TEST_F(contiguous_multimap_test, empty) {
    map_type map{};
    EXPECT_TRUE(std::ranges::empty(map));
    map.insert(value_type{10, 0});
    EXPECT_FALSE(std::ranges::empty(map));
}

TEST_F(contiguous_multimap_test, size) {
    map_type map{};
    EXPECT_EQ(std::ranges::size(map), 0u);
    map.insert(value_type{10, 0});
    EXPECT_EQ(std::ranges::size(map), 1u);
    map.insert(value_type{20, 1});
    EXPECT_EQ(std::ranges::size(map), 2u);
    map.insert(value_type{4, 3});
    EXPECT_EQ(std::ranges::size(map), 3u);
}
