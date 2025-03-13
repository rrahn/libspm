// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

#include <libspm/seqan/alphabet.hpp>

template <typename t>
struct alphabet_test : ::testing::Test
{
    using alphabet_t = t;
};

using test_types = ::testing::Types<
    spm::dna4,
    spm::dna5,
    spm::dna15
>;

TYPED_TEST_SUITE(alphabet_test, test_types);

TYPED_TEST(alphabet_test, alphabet_concept)
{
    EXPECT_TRUE(seqan3::semialphabet<typename TestFixture::alphabet_t>);
    EXPECT_TRUE(seqan3::alphabet<typename TestFixture::alphabet_t>);
    EXPECT_TRUE(seqan3::writable_semialphabet<typename TestFixture::alphabet_t>);
    EXPECT_TRUE(seqan3::writable_alphabet<typename TestFixture::alphabet_t>);
}

TYPED_TEST(alphabet_test, constexpr_alphabet_concept)
{
    EXPECT_TRUE(seqan3::detail::writable_constexpr_semialphabet<typename TestFixture::alphabet_t>);
    EXPECT_TRUE(seqan3::detail::writable_constexpr_alphabet<typename TestFixture::alphabet_t>);
}

TYPED_TEST(alphabet_test, nucleotide_concept)
{
    EXPECT_TRUE(seqan3::nucleotide_alphabet<typename TestFixture::alphabet_t>);
//
}
