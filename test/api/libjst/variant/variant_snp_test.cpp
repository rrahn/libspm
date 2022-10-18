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

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/variant/variant_snp.hpp>

template <typename alphabet_type>
struct snp_test : ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_t = libjst::snp_variant<alphabet_t>;

    snp_t default_snp{};
    snp_t snp{10, seqan3::assign_char_to('C', alphabet_t{})};
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4,
                                    jst::contrib::dna5,
                                    jst::contrib::dna15
                                >;
TYPED_TEST_SUITE(snp_test, test_types);

TYPED_TEST(snp_test, construction)
{
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_t = typename TestFixture::snp_t;

    EXPECT_TRUE((std::is_constructible_v<snp_t, uint32_t, alphabet_t>));
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
    EXPECT_RANGE_EQ(libjst::insertion(this->default_snp), std::vector{seqan3::assign_char_to('A', alphabet_t{})});
    EXPECT_RANGE_EQ(libjst::insertion(this->snp), std::vector{seqan3::assign_char_to('C', alphabet_t{})});
}

TYPED_TEST(snp_test, deletion)
{
    EXPECT_EQ(libjst::deletion(this->default_snp), 1u);
    EXPECT_EQ(libjst::deletion(this->snp), 1u);
}

TYPED_TEST(snp_test, serialise)
{
    using snp_type = typename TestFixture::snp_t;
    using alphabet_t = typename TestFixture::alphabet_t;

    snp_type snp_a_out{0, seqan3::assign_char_to('A', alphabet_t{})};
    snp_type snp_c_out{23, seqan3::assign_char_to('C', alphabet_t{})};
    snp_type snp_g_out{1234, seqan3::assign_char_to('G', alphabet_t{})};
    snp_type snp_t_out{((1 << 30) - 1), seqan3::assign_char_to('T', alphabet_t{})};

    snp_type snp_a_in{};
    snp_type snp_c_in{};
    snp_type snp_g_in{};
    snp_type snp_t_in{};

    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        output_archive(snp_a_out);
        output_archive(snp_c_out);
        output_archive(snp_g_out);
        output_archive(snp_t_out);
    }
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        input_archive(snp_a_in);
        input_archive(snp_c_in);
        input_archive(snp_g_in);
        input_archive(snp_t_in);
    }

    EXPECT_EQ(libjst::position(snp_a_in), libjst::position(snp_a_out));
    EXPECT_RANGE_EQ(libjst::insertion(snp_a_in), libjst::insertion(snp_a_out));

    EXPECT_EQ(libjst::position(snp_c_in), libjst::position(snp_c_out));
    EXPECT_RANGE_EQ(libjst::insertion(snp_c_in), libjst::insertion(snp_c_out));

    EXPECT_EQ(libjst::position(snp_g_in), libjst::position(snp_g_out));
    EXPECT_RANGE_EQ(libjst::insertion(snp_g_in), libjst::insertion(snp_g_out));

    EXPECT_EQ(libjst::position(snp_t_in), libjst::position(snp_t_out));
    EXPECT_RANGE_EQ(libjst::insertion(snp_t_in), libjst::insertion(snp_t_out));
}
