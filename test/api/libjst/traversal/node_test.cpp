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
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_forward.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_model.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>
#include <libjst/sequence_variant/variant_store_covered.hpp>
#include <libjst/traversal/jst_node.hpp>

template <typename alphabet_type>
struct jst_node_test : public ::testing::Test
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
    using fwd_jst_t = libjst::journaled_sequence_tree_forward_<jst_t>;
    using jst_node_t = libjst::jst_node<fwd_jst_t>;

    inline static const std::vector<alphabet_t> base_sequence{seqan3::test::generate_sequence<alphabet_t>(200)};
    inline static const std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

    snp_variant_t snp0{4, seqan3::assign_rank_to(3, alphabet_t{})};
    snp_variant_t snp1{44, seqan3::assign_rank_to(0, alphabet_t{})};
    snp_variant_t snp2{112, seqan3::assign_rank_to(1, alphabet_t{})};
    generic_variant_t var0{44, insertion_sequence, 10};
    generic_variant_t var1{93, insertion_sequence, 0};
    generic_variant_t var2{154, {}, 1};

    jst_t jst;
    fwd_jst_t fwd_jst;

    void SetUp() override
    {
        jst = jst_t{this->base_sequence, 4};

        using value_t = std::ranges::range_value_t<covered_store_t>;

        jst.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}});
        jst.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}});
        jst.insert(value_t{this->snp2, coverage_t{1, 0, 0, 1}});
        jst.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}});
        jst.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}});
        jst.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}});

        fwd_jst = fwd_jst_t{jst};
    }
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(jst_node_test, test_types);

TYPED_TEST(jst_node_test, construction)
{
    using fwd_jst_t = typename TestFixture::fwd_jst_t;
    using jst_node_t = typename TestFixture::jst_node_t;

    EXPECT_TRUE(std::is_default_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_node_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<jst_node_t>);
    EXPECT_TRUE(std::is_destructible_v<jst_node_t>);
    EXPECT_TRUE((std::is_constructible_v<jst_node_t, fwd_jst_t const &, size_t>));
}
