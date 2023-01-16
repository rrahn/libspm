// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <libjst/sequence_tree/node_descriptor.hpp>

struct node_descriptor_test : public ::testing::Test {
};

TEST_F(node_descriptor_test, state_A) {
    libjst::node_descriptor desc{libjst::node_state::branching_after_left_end};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::branching_after_left_end);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_TRUE(desc.is_branching());

    EXPECT_FALSE(desc.left_break().from_left_begin());
    EXPECT_TRUE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_TRUE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_B) {
    libjst::node_descriptor desc{libjst::node_state::last_branching_after_left_end};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::last_branching_after_left_end);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_TRUE(desc.is_branching());

    EXPECT_FALSE(desc.left_break().from_left_begin());
    EXPECT_TRUE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_TRUE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_E) {
    libjst::node_descriptor desc{libjst::node_state::branching_after_left_begin};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::branching_after_left_begin);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_TRUE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_TRUE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_F) {
    libjst::node_descriptor desc{libjst::node_state::last_branching_after_left_begin};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::last_branching_after_left_begin);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_TRUE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_TRUE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_C) {
    libjst::node_descriptor desc{libjst::node_state::last_non_branching_left_only};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::last_non_branching_left_only);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_FALSE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_TRUE(desc.right_break().from_left_end());
    EXPECT_FALSE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_D) {
    libjst::node_descriptor desc{libjst::node_state::non_branching_left_only};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::non_branching_left_only);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_FALSE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_TRUE(desc.right_break().from_left_end());
    EXPECT_FALSE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_G) {
    libjst::node_descriptor desc{libjst::node_state::non_branching_after_left};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::non_branching_after_left);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_FALSE(desc.is_branching());

    EXPECT_FALSE(desc.left_break().from_left_begin());
    EXPECT_TRUE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_FALSE(desc.right_break().from_right_begin());
    EXPECT_TRUE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_H) {
    libjst::node_descriptor desc{libjst::node_state::non_branching_including_left};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::non_branching_including_left);

    EXPECT_TRUE(desc.from_reference());
    EXPECT_FALSE(desc.from_variant());
    EXPECT_FALSE(desc.on_alternate_path());
    EXPECT_FALSE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_FALSE(desc.right_break().from_left_end());
    EXPECT_FALSE(desc.right_break().from_right_begin());
    EXPECT_TRUE(desc.right_break().from_right_end());
}

TEST_F(node_descriptor_test, state_variant) {
    libjst::node_descriptor desc{libjst::node_state::variant};

    EXPECT_EQ(static_cast<libjst::node_state>(desc), libjst::node_state::variant);

    EXPECT_FALSE(desc.from_reference());
    EXPECT_TRUE(desc.from_variant());
    EXPECT_TRUE(desc.on_alternate_path());
    EXPECT_FALSE(desc.is_branching());

    EXPECT_TRUE(desc.left_break().from_left_begin());
    EXPECT_FALSE(desc.left_break().from_left_end());
    EXPECT_FALSE(desc.left_break().from_right_begin());
    EXPECT_FALSE(desc.left_break().from_right_end());

    EXPECT_FALSE(desc.right_break().from_left_begin());
    EXPECT_TRUE(desc.right_break().from_left_end());
    EXPECT_FALSE(desc.right_break().from_right_begin());
    EXPECT_FALSE(desc.right_break().from_right_end());
}
