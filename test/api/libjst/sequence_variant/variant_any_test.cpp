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

#include <libjst/sequence_variant/variant_any.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>

template <typename alphabet_type>
struct any_variant_test : ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;

    using position_t = std::common_type_t<libjst::variant_position_t<snp_variant_t>, libjst::variant_position_t<generic_variant_t>>;
    using insertion_t = std::common_type_t<libjst::variant_insertion_t<snp_variant_t>, libjst::variant_insertion_t<generic_variant_t>>;
    using deletion_t = std::common_type_t<libjst::variant_deletion_t<snp_variant_t>, libjst::variant_deletion_t<generic_variant_t>>;
    using any_variant_t = libjst::any_variant<position_t, insertion_t, deletion_t>;

    inline static const std::vector<alphabet_t> insertion_sequence = seqan3::test::generate_sequence<alphabet_t>(10);

    snp_variant_t snp_var{6, seqan3::assign_rank_to(2, alphabet_t{})};
    generic_variant_t variant_sub{10, insertion_sequence, 10};
    generic_variant_t variant_ins{20, insertion_sequence, 0};
    generic_variant_t variant_del{34, {}, 7};
};

using test_types = ::testing::Types<//jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(any_variant_test, test_types);

TYPED_TEST(any_variant_test, construction)
{
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;
    using any_variant_t = typename TestFixture::any_variant_t;

    EXPECT_TRUE((std::is_constructible_v<any_variant_t, snp_variant_t>));
    EXPECT_TRUE((std::is_constructible_v<any_variant_t, generic_variant_t>));
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<any_variant_t>);
    EXPECT_FALSE(std::is_copy_constructible_v<any_variant_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<any_variant_t>);
    EXPECT_FALSE(std::is_copy_assignable_v<any_variant_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<any_variant_t>);
    EXPECT_TRUE(std::is_destructible_v<any_variant_t>);
}

TYPED_TEST(any_variant_test, concept)
{
    using any_variant_t = typename TestFixture::any_variant_t;
    EXPECT_TRUE(libjst::sequence_variant<any_variant_t>);
    EXPECT_TRUE(libjst::sequence_variant<any_variant_t &>);
    EXPECT_TRUE(libjst::sequence_variant<any_variant_t const>);
    EXPECT_TRUE(libjst::sequence_variant<any_variant_t const &>);
}

TYPED_TEST(any_variant_test, position)
{
    using any_variant_t = typename TestFixture::any_variant_t;
    EXPECT_EQ(libjst::position(any_variant_t{this->snp_var}), libjst::position(this->snp_var));
    EXPECT_EQ(libjst::position(any_variant_t{this->variant_sub}), libjst::position(this->variant_sub));
    EXPECT_EQ(libjst::position(any_variant_t{this->variant_ins}), libjst::position(this->variant_ins));
    EXPECT_EQ(libjst::position(any_variant_t{this->variant_del}), libjst::position(this->variant_del));
}

TYPED_TEST(any_variant_test, insertion)
{
    using any_variant_t = typename TestFixture::any_variant_t;
    EXPECT_RANGE_EQ(libjst::insertion(any_variant_t{this->snp_var}), libjst::insertion(this->snp_var));
    EXPECT_RANGE_EQ(libjst::insertion(any_variant_t{this->variant_sub}), libjst::insertion(this->variant_sub));
    EXPECT_RANGE_EQ(libjst::insertion(any_variant_t{this->variant_ins}), libjst::insertion(this->variant_ins));
    EXPECT_RANGE_EQ(libjst::insertion(any_variant_t{this->variant_del}), libjst::insertion(this->variant_del));
}

TYPED_TEST(any_variant_test, deletion)
{
    using any_variant_t = typename TestFixture::any_variant_t;
    EXPECT_EQ(libjst::deletion(any_variant_t{this->snp_var}), libjst::deletion(this->snp_var));
    EXPECT_EQ(libjst::deletion(any_variant_t{this->variant_sub}), libjst::deletion(this->variant_sub));
    EXPECT_EQ(libjst::deletion(any_variant_t{this->variant_ins}), libjst::deletion(this->variant_ins));
    EXPECT_EQ(libjst::deletion(any_variant_t{this->variant_del}), libjst::deletion(this->variant_del));
}
