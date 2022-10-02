// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/sequence_variant/variant_snp.hpp>

template <typename alphabet_type>
struct snp_test : ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_t = libjst::snp_variant<alphabet_t>;

    snp_t default_snp{};
    snp_t snp{10, seqan3::assign_rank_to(1, alphabet_t{})};
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(snp_test, test_types);

TYPED_TEST(snp_test, construction)
{
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_t = typename TestFixture::snp_t;

    EXPECT_TRUE((std::is_nothrow_constructible_v<snp_t, uint32_t, alphabet_t>));
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<snp_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<snp_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<snp_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<snp_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<snp_t>);
    EXPECT_TRUE(std::is_destructible_v<snp_t>);
}

TYPED_TEST(snp_test, concept)
{
    using snp_t = typename TestFixture::snp_t;
    EXPECT_TRUE(libjst::sequence_variant<snp_t>);
    EXPECT_TRUE(libjst::sequence_variant<snp_t &>);
    EXPECT_TRUE(libjst::sequence_variant<snp_t const>);
    EXPECT_TRUE(libjst::sequence_variant<snp_t const &>);
}

TYPED_TEST(snp_test, sizeof)
{
    using snp_t = typename TestFixture::snp_t;
    EXPECT_EQ(sizeof(snp_t), 4);
}

TYPED_TEST(snp_test, position)
{
    EXPECT_EQ(libjst::position(this->default_snp), 0u);
    EXPECT_EQ(libjst::position(this->snp), 10u);
}

TYPED_TEST(snp_test, insertion)
{
    using alphabet_t = typename TestFixture::alphabet_t;
    EXPECT_RANGE_EQ(libjst::insertion(this->default_snp), std::vector{seqan3::assign_rank_to(0, alphabet_t{})});
    EXPECT_RANGE_EQ(libjst::insertion(this->snp), std::vector{seqan3::assign_rank_to(1, alphabet_t{})});
}

TYPED_TEST(snp_test, deletion)
{
    EXPECT_EQ(libjst::deletion(this->default_snp), 1u);
    EXPECT_EQ(libjst::deletion(this->snp), 1u);
}
