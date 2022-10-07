// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/sequence_variant/variant_generic.hpp>

template <typename alphabet_type>
struct generic_variant_test : ::testing::Test
{

    using alphabet_t = alphabet_type;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;

    inline static const std::vector<alphabet_t> insertion_sequence = seqan3::test::generate_sequence<alphabet_t>(10);

    generic_variant_t default_var{};
    generic_variant_t variant_sub{10, insertion_sequence, 10};
    generic_variant_t variant_ins{20, insertion_sequence, 0};
    generic_variant_t variant_del{34, {}, 7};
};

using test_types = ::testing::Types<//jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(generic_variant_test, test_types);

TYPED_TEST(generic_variant_test, construction)
{
    using alphabet_t = typename TestFixture::alphabet_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    EXPECT_TRUE((std::is_nothrow_constructible_v<generic_variant_t, uint32_t, std::vector<alphabet_t>, uint32_t>));
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<generic_variant_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<generic_variant_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<generic_variant_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<generic_variant_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<generic_variant_t>);
    EXPECT_TRUE(std::is_destructible_v<generic_variant_t>);
}

TYPED_TEST(generic_variant_test, concept)
{
    using generic_variant_t = typename TestFixture::generic_variant_t;
    EXPECT_TRUE(libjst::sequence_variant<generic_variant_t>);
    EXPECT_TRUE(libjst::sequence_variant<generic_variant_t &>);
    EXPECT_TRUE(libjst::sequence_variant<generic_variant_t const>);
    EXPECT_TRUE(libjst::sequence_variant<generic_variant_t const &>);
}

TYPED_TEST(generic_variant_test, position)
{
    EXPECT_EQ(libjst::position(this->default_var), 0u);
    EXPECT_EQ(libjst::position(this->variant_sub), 10u);
    EXPECT_EQ(libjst::position(this->variant_ins), 20u);
    EXPECT_EQ(libjst::position(this->variant_del), 34u);
}

TYPED_TEST(generic_variant_test, insertion)
{
    using alphabet_t = typename TestFixture::alphabet_t;
    EXPECT_RANGE_EQ(libjst::insertion(this->default_var), std::vector<alphabet_t>{});
    EXPECT_RANGE_EQ(libjst::insertion(this->variant_sub), TestFixture::insertion_sequence);
    EXPECT_RANGE_EQ(libjst::insertion(this->variant_ins), TestFixture::insertion_sequence);
    EXPECT_RANGE_EQ(libjst::insertion(this->variant_del), std::vector<alphabet_t>{});
}

TYPED_TEST(generic_variant_test, deletion)
{
    EXPECT_EQ(libjst::deletion(this->default_var), 0u);
    EXPECT_EQ(libjst::deletion(this->variant_sub), std::ranges::size(TestFixture::insertion_sequence));
    EXPECT_EQ(libjst::deletion(this->variant_ins), 0u);
    EXPECT_EQ(libjst::deletion(this->variant_del), 7u);
}

TYPED_TEST(generic_variant_test, serialise)
{
    using generic_variant_t = typename TestFixture::generic_variant_t;
    using alphabet_t = typename TestFixture::alphabet_t;

    generic_variant_t var_sub_out{0, this->insertion_sequence, (uint32_t)std::ranges::size(this->insertion_sequence)};
    generic_variant_t var_del_out{1234, std::vector<alphabet_t>{}, 15};
    generic_variant_t var_ins_out{((1 << 30) - 1), this->insertion_sequence, 0};

    generic_variant_t var_sub_in{};
    generic_variant_t var_ins_in{};
    generic_variant_t var_del_in{};

    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        output_archive(var_sub_out);
        output_archive(var_del_out);
        output_archive(var_ins_out);
    }
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        input_archive(var_sub_in);
        input_archive(var_del_in);
        input_archive(var_ins_in);
    }

    EXPECT_EQ(libjst::position(var_sub_in), libjst::position(var_sub_out));
    EXPECT_EQ(libjst::deletion(var_sub_in), libjst::deletion(var_sub_out));
    EXPECT_RANGE_EQ(libjst::insertion(var_sub_in), libjst::insertion(var_sub_out));

    EXPECT_EQ(libjst::position(var_del_in), libjst::position(var_del_out));
    EXPECT_EQ(libjst::deletion(var_del_in), libjst::deletion(var_del_out));
    EXPECT_RANGE_EQ(libjst::insertion(var_del_in), libjst::insertion(var_del_out));

    EXPECT_EQ(libjst::position(var_ins_in), libjst::position(var_ins_out));
    EXPECT_EQ(libjst::deletion(var_ins_in), libjst::deletion(var_ins_out));
    EXPECT_RANGE_EQ(libjst::insertion(var_ins_in), libjst::insertion(var_ins_out));
}
