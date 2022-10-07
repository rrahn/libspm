// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/journaled_sequence_tree/concept.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_model.hpp>
#include <libjst/journaled_sequence_tree/journaled_sequence_tree_forward.hpp>
#include <libjst/utility/bit_vector.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>
#include <libjst/sequence_variant/variant_store_covered.hpp>

#include <libjst/journaled_sequence_tree/serialiser_direct.hpp>
#include <libjst/journaled_sequence_tree/serialiser_delegate.hpp>

TEST(journaled_sequence_tree_serialiser_test, protoype_jst)
{
    using alphabet_t = jst::contrib::dna4;
    using sequence_t = std::vector<alphabet_t>;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;
    using coverage_t = libjst::bit_vector<>;

    using snp_store_t = std::vector<snp_variant_t>;
    using generic_store_t = std::vector<generic_variant_t>;
    using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;
    using covered_store_t = libjst::variant_store_covered<composite_store_t, libjst::bit_vector<>>;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    using jst_t = libjst::journaled_sequence_tree_model<sequence_t, covered_store_t>;
    using fwd_jst_t = libjst::journaled_sequence_tree_forward_<jst_t>;

    std::vector<alphabet_t> base_sequence{seqan3::test::generate_sequence<alphabet_t>(200)};
    std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

    snp_variant_t snp0{4, seqan3::assign_rank_to(3, alphabet_t{})};
    snp_variant_t snp1{112, seqan3::assign_rank_to(0, alphabet_t{})};
    generic_variant_t var0{44, insertion_sequence, 10};
    generic_variant_t var1{93, insertion_sequence, 0};
    generic_variant_t var2{154, {}, 1};

    jst_t jst_out{base_sequence, 4};

    EXPECT_TRUE((jst_out.insert(value_t{snp0, coverage_t{0, 0, 0, 1}})));
    EXPECT_TRUE((jst_out.insert(value_t{var0, coverage_t{0, 0, 1, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{var1, coverage_t{0, 1, 0, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{snp1, coverage_t{1, 0, 0, 0}})));
    EXPECT_TRUE((jst_out.insert(value_t{var2, coverage_t{0, 0, 1, 1}})));

    // test output stream and input stream in same buffer
    fwd_jst_t fwd_jst_out{jst_out};
    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        auto arch = output_archive | libjst::direct_serialiser(base_sequence)
                                   | libjst::delegate_serialiser(jst_out);
        libjst::save(fwd_jst_out, arch);
    }

    std::vector<alphabet_t> base_sequence_in{};
    jst_t jst_in{base_sequence_in, 0}; // some default values.
    fwd_jst_t fwd_jst_in{jst_in};
    { // input arcihve
        cereal::JSONInputArchive input_archive(archive_stream);
        auto arch = input_archive | libjst::direct_serialiser(base_sequence_in)
                                  | libjst::delegate_serialiser(jst_in);
        libjst::load(fwd_jst_in, arch);
    }

    EXPECT_RANGE_EQ(libjst::base_sequence(fwd_jst_in), libjst::base_sequence(fwd_jst_out));
    EXPECT_EQ(libjst::size(jst_in), libjst::size(jst_out));
    auto const & variant_store_out = libjst::variant_store(fwd_jst_out);
    auto const & variant_store_in = libjst::variant_store(fwd_jst_in);
    EXPECT_EQ(std::ranges::size(variant_store_in), std::ranges::size(variant_store_out));
    for (unsigned i = 0; i < std::ranges::size(variant_store_in); ++i)
    {
        EXPECT_EQ(libjst::position(variant_store_in[i]), libjst::position(variant_store_out[i]));
        EXPECT_EQ(libjst::deletion(variant_store_in[i]), libjst::deletion(variant_store_out[i]));
        EXPECT_RANGE_EQ(libjst::insertion(variant_store_in[i]), libjst::insertion(variant_store_out[i]));
        EXPECT_RANGE_EQ(libjst::coverage(variant_store_in[i]), libjst::coverage(variant_store_out[i]));
    }
}
