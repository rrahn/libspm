// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <jstmap/search/seed_prefix_seek_position.hpp>

struct seed_prefix_seek_position_test : public ::testing::Test
{
        //   * *   *     *
      // 0 1 2 3 4 5 6 7 8 9
      // 9 8 7 6 5 4 3 2 1 0
        //0123456789012345
    size_t breakends_count = 10;
};

TEST_F(seed_prefix_seek_position_test, reference_position) {
    {
        libjst::seek_position seed_position{};
        seed_position.reset(4, libjst::breakpoint_end::low);
        jstmap::seed_prefix_seek_position position{std::move(seed_position), breakends_count};

        EXPECT_EQ(position.get_variant_index(), 4);
        position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
            if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
                EXPECT_EQ(descriptor, libjst::breakpoint_end::low);
            } else {
                FAIL() << "Expected node from reference path.";
            }
        });
    }

    {
        libjst::seek_position seed_position{};
        seed_position.reset(2, libjst::breakpoint_end::high);
        jstmap::seed_prefix_seek_position position{std::move(seed_position), breakends_count};

        EXPECT_EQ(position.get_variant_index(), 6);
        position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
            if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
                EXPECT_EQ(descriptor, libjst::breakpoint_end::high);
            } else {
                FAIL() << "Expected node from reference path.";
            }
        });
    }
}

TEST_F(seed_prefix_seek_position_test, alternate_position) {
    {
        libjst::seek_position seed_position{};
        seed_position.initiate_alternate_node(4);
        seed_position.next_alternate_node(true);
        seed_position.next_alternate_node(false);
        seed_position.next_alternate_node(false);
        seed_position.next_alternate_node(true);
        jstmap::seed_prefix_seek_position position{std::move(seed_position), breakends_count};

        EXPECT_EQ(position.get_variant_index(), 5);
        position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
            if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
                FAIL() << "Expected node from alternate path.";
            } else {
                EXPECT_EQ(descriptor.size(), 1u);
            }
        });
    }

    {
        libjst::seek_position seed_position{};
        seed_position.initiate_alternate_node(2);
        seed_position.next_alternate_node(true);
        seed_position.next_alternate_node(true);
        seed_position.next_alternate_node(false);
        seed_position.next_alternate_node(false);

        jstmap::seed_prefix_seek_position position{std::move(seed_position), breakends_count};

        EXPECT_EQ(position.get_variant_index(), 7);
        position.visit([&] <typename descriptor_t> (descriptor_t descriptor) {
            if constexpr (std::same_as<descriptor_t, libjst::breakpoint_end>) {
                FAIL() << "Expected node from alternate path.";
            } else {
                EXPECT_EQ(descriptor.size(), 1u);
            }
        });
    }
}
