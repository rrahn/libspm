// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>
#include <span>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/utility/bit_vector.hpp>

#include <libjst/journal.hpp>
#include <libjst/variant/variant_snp.hpp>
#include <libjst/variant/variant_store_covered.hpp>
#include <libjst/tree/jst_node_base.hpp>
#include <libjst/tree/branch_state.hpp>

template <typename alphabet_type>
struct jst_node_test : public ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_t = libjst::snp_variant<alphabet_t>;
    using position_t = libjst::variant_position_t<snp_t>;
    using coverage_t = libjst::bit_vector<>;
    using snp_store_t = std::vector<snp_t>;
    using store_t = libjst::variant_store_covered<snp_store_t, coverage_t>;

    using journal_t = libjst::journal<position_t, std::span<alphabet_t>>;
    using node_value_t = libjst::jst_node_value<journal_t, coverage_t>;
    using store_iterator_t = std::ranges::iterator_t<store_t const>;
    using jst_node_t = libjst::jst_node_base<node_value_t, store_iterator_t>;
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    jst::contrib::dna5
                                    >;
TYPED_TEST_SUITE(jst_node_test, test_types);

TYPED_TEST(jst_node_test, construction)
{
    using node_value_t = typename TestFixture::node_value_t;
    using store_iterator_t = typename TestFixture::store_iterator_t;
    using jst_node_t = typename TestFixture::jst_node_t;

    EXPECT_TRUE(std::is_default_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<jst_node_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_node_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<jst_node_t>);
    EXPECT_TRUE(std::is_destructible_v<jst_node_t>);
    EXPECT_TRUE((std::is_constructible_v<jst_node_t, node_value_t, store_iterator_t, store_iterator_t, size_t>));
}
