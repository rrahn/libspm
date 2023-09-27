// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <algorithm>

#include <libjst/sequence_tree/journaled_sequence_label.hpp>

#include "rcs_store_mock.hpp"

struct js_label_fixture : public ::testing::Test {

    using source_t = std::string_view;
    using source_view_t = libjst::detail::subrange_t<source_t const &>;
    using label_type = libjst::journaled_sequence_label<int32_t, source_view_t>;
                                    //0         1
                                    //01234567890123
    static constexpr source_t source{"garfieldthecat"};

};

TEST_F(js_label_fixture, create_from_source) {
    label_type lbl{this->source};
    EXPECT_EQ(lbl.get_left_position(), 0);
    EXPECT_EQ(lbl.get_right_position(), std::ranges::ssize(this->source));
    EXPECT_EQ(lbl.label_size(), std::ranges::ssize(this->source));
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), this->source));
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), this->source));
}

TEST_F(js_label_fixture, record_variant) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.record_variant(jst::test::variant{.position = 8, .insertion = ""sv, .deletion = 3});
    EXPECT_EQ(lbl.get_left_position(), 8);
    EXPECT_EQ(lbl.get_right_position(), 8);
    EXPECT_EQ(lbl.label_size(), 0u);
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldcat"sv));
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), ""sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "fat"sv, .deletion = 0});
    EXPECT_EQ(lbl.get_left_position(), 8);
    EXPECT_EQ(lbl.get_right_position(), 11);
    EXPECT_EQ(lbl.label_size(), 3u);
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldfatcat"sv));
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "fat"sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "r"sv, .deletion = 1});
    EXPECT_EQ(lbl.get_left_position(), 11);
    EXPECT_EQ(lbl.get_right_position(), 12);
    EXPECT_EQ(lbl.label_size(), 1u);
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldfatrat"sv));
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "r"sv));
}

TEST_F(js_label_fixture, path_sequence) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.record_variant(jst::test::variant{.position = 8, .insertion = ""sv, .deletion = 3});
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldcat"sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "fat"sv, .deletion = 0});
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldfatcat"sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "r"sv, .deletion = 1});
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), "garfieldfatrat"sv));
}

TEST_F(js_label_fixture, node_sequence) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.record_variant(jst::test::variant{.position = 8, .insertion = ""sv, .deletion = 3});
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), ""sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "fat"sv, .deletion = 0});
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "fat"sv));

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "r"sv, .deletion = 1});
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "r"sv));
}

TEST_F(js_label_fixture, get_left_position) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.record_variant(jst::test::variant{.position = 8, .insertion = ""sv, .deletion = 3});
    EXPECT_EQ(lbl.get_left_position(), 8);

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "fat"sv, .deletion = 0});
    EXPECT_EQ(lbl.get_left_position(), 8);

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "r"sv, .deletion = 1});
    EXPECT_EQ(lbl.get_left_position(), 11);
}

TEST_F(js_label_fixture, get_right_position) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.record_variant(jst::test::variant{.position = 8, .insertion = ""sv, .deletion = 3});
    EXPECT_EQ(lbl.get_right_position(), 8);

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "fat"sv, .deletion = 0});
    EXPECT_EQ(lbl.get_right_position(), 11);

    lbl.record_variant(jst::test::variant{.position = 11, .insertion = "r"sv, .deletion = 1});
    EXPECT_EQ(lbl.get_right_position(), 12);
}

TEST_F(js_label_fixture, reset_positions) {
    using namespace std::literals;
    label_type lbl{this->source};

    lbl.reset_positions(0, 0);
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), ""sv));
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), this->source));

    lbl.reset_positions(8, 11);
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "the"sv));
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), this->source));

    lbl.reset_positions(11, 14);
    EXPECT_TRUE(std::ranges::equal(lbl.node_sequence(), "cat"sv));
    EXPECT_TRUE(std::ranges::equal(lbl.path_sequence(), this->source));
}
