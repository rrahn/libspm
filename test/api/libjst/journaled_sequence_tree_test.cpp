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

#include <cereal/archives/json.hpp> // output archive for testing
#include <cereal/types/string.hpp> // sereialise std::string

#include <seqan3/alphabet/adaptation/char.hpp> // allow std::string be recognised as seqan3::sequence

#include <libjst/journaled_sequence_tree.hpp>

struct journaled_sequence_tree_fixture : public ::testing::Test
{
    using sequence_t = std::string;
    using jst_t = libjst::journaled_sequence_tree<sequence_t>;
    using alignment_t = std::pair<sequence_t, sequence_t>;

    sequence_t reference{"aaaabbbbcccc"};

    alignment_t alignment1{"aaaabbbbcccc-----", "------------aabbcc"};
    alignment_t alignment2{"aaaabbbbcccc-----", "------------abcabc"};
    alignment_t alignment3{"aaaa--bbbb--cccc--", "----cc----aa----bb"};
};

TEST_F(journaled_sequence_tree_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_move_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_move_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<jst_t>);
    EXPECT_TRUE((std::is_constructible_v<jst_t, sequence_t>)); // Set explicit reference sequence.
    EXPECT_TRUE((std::is_constructible_v<jst_t, sequence_t &&>)); // Set explicit reference sequence.
    // I might be able to convert a journaled sequence to an alignment.
    // TODO: Construct from multiple sequence alignment: reference is consensus sequence
}

TEST_F(journaled_sequence_tree_fixture, reference)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    EXPECT_EQ(jst.reference(), reference);
}

TEST_F(journaled_sequence_tree_fixture, size)
{
    jst_t jst{std::move(reference)};

    EXPECT_EQ(jst.size(), 0u);
}

TEST_F(journaled_sequence_tree_fixture, add)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 1u);

    jst.add(alignment2); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 2u);

    jst.add(alignment3); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 3u);

    alignment_t alignment_wrong_reference{"aaaabbbbccc-----x", alignment1.second};
    EXPECT_THROW(jst.add(alignment_wrong_reference), std::invalid_argument);

    alignment_t alignment_wrong_order{alignment1.second, alignment1.first};
    EXPECT_THROW(jst.add(alignment_wrong_order), std::invalid_argument);
}

TEST_F(journaled_sequence_tree_fixture, save)
{

    std::stringstream output_stream{};
    cereal::JSONOutputArchive output_archive(output_stream);

    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1);
    jst.add(alignment2);
    jst.add(alignment3);

    jst.save(output_archive);

    std::string_view expected_output = R"json({
    "value0": "aaaabbbbcccc",
    "value1": [
        "aabbcc",
        "abcabc",
        "ccaabb"
    ])json";

    EXPECT_EQ(output_stream.str(), expected_output);
}
