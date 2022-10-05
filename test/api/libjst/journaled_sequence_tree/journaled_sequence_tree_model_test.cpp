// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_model.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>
#include <libjst/sequence_variant/variant_store_covered.hpp>

template <typename alphabet_type>
struct journaled_sequence_tree_model_test : public ::testing::Test
{
    using alphabet_t = alphabet_type;
    using sequence_t = std::vector<alphabet_t>;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;
    using coverage_t = libjst::bit_vector<>;

    using snp_store_t = std::vector<snp_variant_t>;
    using generic_store_t = std::vector<generic_variant_t>;
    using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;
    using covered_store_t = libjst::variant_store_covered<composite_store_t, libjst::bit_vector<>>;

    using jst_t = libjst::journaled_sequence_tree_model<sequence_t, covered_store_t>;

    inline static const std::vector<alphabet_t> base_sequence{seqan3::test::generate_sequence<alphabet_t>(200)};
    inline static const std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

    snp_variant_t snp0{4, seqan3::assign_rank_to(3, alphabet_t{})};
    snp_variant_t snp1{112, seqan3::assign_rank_to(0, alphabet_t{})};
    generic_variant_t var0{44, insertion_sequence, 10};
    generic_variant_t var1{93, insertion_sequence, 0};
    generic_variant_t var2{154, {}, 1};
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(journaled_sequence_tree_model_test, test_types);

TYPED_TEST(journaled_sequence_tree_model_test, construction)
{
    using jst_t = typename TestFixture::jst_t;
    using sequence_t = typename TestFixture::sequence_t;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_destructible_v<jst_t>);
    EXPECT_TRUE((std::is_constructible_v<jst_t, sequence_t, size_t>));
}

TYPED_TEST(journaled_sequence_tree_model_test, concept)
{
    using jst_t = typename TestFixture::jst_t;

    EXPECT_TRUE(libjst::journaled_sequence_tree<jst_t>);
    EXPECT_TRUE(libjst::journaled_sequence_tree<jst_t &>);
    EXPECT_TRUE(libjst::journaled_sequence_tree<jst_t const &>);
    EXPECT_FALSE(libjst::traversable_journaled_sequence_tree<jst_t>);
    EXPECT_FALSE(libjst::traversable_journaled_sequence_tree<jst_t &>);
    EXPECT_FALSE(libjst::traversable_journaled_sequence_tree<jst_t const &>);
}

TYPED_TEST(journaled_sequence_tree_model_test, insert)
{
    using jst_t = typename TestFixture::jst_t;
    using covered_store_t = typename TestFixture::covered_store_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>; // does not have this in the jst!

    jst_t jst{this->base_sequence, 4};

    EXPECT_TRUE((jst.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}})));
    EXPECT_TRUE((jst.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}})));
    EXPECT_TRUE((jst.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}})));
    EXPECT_TRUE((jst.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}})));
    EXPECT_TRUE((jst.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}})));
}

TYPED_TEST(journaled_sequence_tree_model_test, size)
{
    using jst_t = typename TestFixture::jst_t;

    jst_t def_jst{};
    jst_t jst{this->base_sequence, 4};

    EXPECT_EQ(libjst::size(def_jst), 0);
    EXPECT_EQ(libjst::size(jst), 4u);
}
