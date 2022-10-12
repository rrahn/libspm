// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/utility/bit_vector.hpp>
#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>
#include <libjst/sequence_variant/variant_store_covered.hpp>

template <typename alphabet_type>
struct variant_store_covered_test : public ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;
    using coverage_t = libjst::bit_vector<>;

    using snp_store_t = std::vector<snp_variant_t>;
    using generic_store_t = std::vector<generic_variant_t>;
    using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;
    using covered_store_t = libjst::variant_store_covered<composite_store_t, libjst::bit_vector<>>;

    inline static const std::vector<alphabet_t> insertion_sequence{seqan3::test::generate_sequence<alphabet_t>(10)};

    snp_variant_t snp0{4, seqan3::assign_rank_to(3, alphabet_t{})};
    snp_variant_t snp1{112, seqan3::assign_rank_to(0, alphabet_t{})};
    generic_variant_t var0{44, insertion_sequence, 10};
    generic_variant_t var1{93, insertion_sequence, 0};
    generic_variant_t var2{154, {}, 1};

    auto make_insertion()
    {
        return seqan3::test::generate_sequence<alphabet_t>(10);
    }
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4
                                    >;
TYPED_TEST_SUITE(variant_store_covered_test, test_types);

TYPED_TEST(variant_store_covered_test, construction)
{
    using covered_store_t = typename TestFixture::covered_store_t;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<covered_store_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<covered_store_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<covered_store_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<covered_store_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<covered_store_t>);
    EXPECT_TRUE(std::is_destructible_v<covered_store_t>);
}

TYPED_TEST(variant_store_covered_test, concept)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    EXPECT_TRUE(std::ranges::random_access_range<covered_store_t>);
    EXPECT_TRUE(libjst::sequence_variant_store<covered_store_t>);
    EXPECT_TRUE(libjst::sequence_variant_store<covered_store_t &>);
    EXPECT_TRUE(libjst::sequence_variant_store<covered_store_t const &>);
    EXPECT_TRUE(libjst::covered_sequence_variant_store<covered_store_t>);
    EXPECT_TRUE(libjst::covered_sequence_variant_store<covered_store_t &>);
    EXPECT_TRUE(libjst::covered_sequence_variant_store<covered_store_t const &>);
}

TYPED_TEST(variant_store_covered_test, type_traits)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using coverage_t = typename TestFixture::coverage_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;
    using reference_t = std::ranges::range_reference_t<covered_store_t>;
    using const_reference_t = std::ranges::range_reference_t<covered_store_t const>;

    EXPECT_TRUE(libjst::covered_sequence_variant<value_t>);
    EXPECT_TRUE(libjst::covered_sequence_variant<reference_t>);
    // EXPECT_TRUE((std::constructible_from<reference_t, value_t &>));
    EXPECT_TRUE((std::constructible_from<const_reference_t, value_t const &>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t>));

    snp_variant_t snp{4, seqan3::assign_rank_to(3, alphabet_t{})};
    coverage_t coverage{0, 1, 1, 0};
    reference_t ref{snp, coverage};
    value_t val{snp, coverage};
    EXPECT_EQ(libjst::position(snp), 4u);
    EXPECT_EQ(libjst::position(ref), 4u);
    EXPECT_EQ(libjst::position(val), 4u);
    EXPECT_RANGE_EQ(libjst::coverage(ref), (coverage_t{0, 1, 1, 0}));
    EXPECT_RANGE_EQ(libjst::coverage(val), (coverage_t{0, 1, 1, 0}));

    snp = snp_variant_t{10, seqan3::assign_rank_to(0, alphabet_t{})};
    coverage[0] = 1;
    coverage[1] = 0;
    EXPECT_EQ(libjst::position(snp), 10u);
    EXPECT_EQ(libjst::position(ref), 10u);
    EXPECT_EQ(libjst::position(val), 4u);
    EXPECT_RANGE_EQ(libjst::coverage(ref), (coverage_t{1, 0, 1, 0}));
    EXPECT_RANGE_EQ(libjst::coverage(val), (coverage_t{0, 1, 1, 0}));
}

TYPED_TEST(variant_store_covered_test, insert)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    covered_store_t store{};

    EXPECT_EQ((store.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}}) - store.begin()), 0);
    EXPECT_EQ((store.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}}) - store.begin()), 2);
    EXPECT_EQ((store.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}}) - store.begin()), 4);
}

TYPED_TEST(variant_store_covered_test, emplace)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using coverage_t = typename TestFixture::coverage_t;

    covered_store_t store{};

    EXPECT_EQ((store.emplace(this->snp0, coverage_t{0, 0, 0, 1}) - store.begin()), 0);
    EXPECT_EQ((store.emplace(this->var0, coverage_t{0, 0, 1, 0}) - store.begin()), 1);
    EXPECT_EQ((store.emplace(this->var1, coverage_t{0, 1, 0, 0}) - store.begin()), 2);
    EXPECT_EQ((store.emplace(this->snp1, coverage_t{1, 0, 0, 0}) - store.begin()), 1);
    EXPECT_EQ((store.emplace(this->var2, coverage_t{0, 0, 1, 1}) - store.begin()), 4);
}

TYPED_TEST(variant_store_covered_test, size)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    covered_store_t store{};
    coverage_t cov{0, 1, 0, 1};
    EXPECT_EQ(store.size(), 0u);
    EXPECT_NO_THROW((store.insert(value_t{this->snp0, cov})));
    EXPECT_EQ(store.size(), 1u);
    EXPECT_NO_THROW((store.insert(value_t{this->var0, cov})));
    EXPECT_EQ(store.size(), 2u);
    EXPECT_NO_THROW((store.insert(value_t{this->var1, cov})));
    EXPECT_EQ(store.size(), 3u);
    EXPECT_NO_THROW((store.insert(value_t{this->snp1, cov})));
    EXPECT_EQ(store.size(), 4u);
    EXPECT_NO_THROW((store.insert(value_t{this->var2, cov})));
    EXPECT_EQ(store.size(), 5u);
}

TYPED_TEST(variant_store_covered_test, subscript)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    covered_store_t store{};

    EXPECT_EQ((store.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}}) - store.begin()), 0);
    EXPECT_EQ((store.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}}) - store.begin()), 2);
    EXPECT_EQ((store.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}}) - store.begin()), 4);

    EXPECT_EQ(libjst::position(store[0]), 4u);
    EXPECT_EQ(libjst::position(store[1]), 112u);
    EXPECT_EQ(libjst::position(store[2]), 44u);
    EXPECT_EQ(libjst::position(store[3]), 93u);
    EXPECT_EQ(libjst::position(store[4]), 154u);

    EXPECT_RANGE_EQ(libjst::insertion(store[0]), (std::vector{seqan3::assign_rank_to(3, alphabet_t{})}));
    EXPECT_RANGE_EQ(libjst::insertion(store[1]), (std::vector{seqan3::assign_rank_to(0, alphabet_t{})}));
    EXPECT_RANGE_EQ(libjst::insertion(store[2]), this->insertion_sequence);
    EXPECT_RANGE_EQ(libjst::insertion(store[3]), this->insertion_sequence);
    EXPECT_RANGE_EQ(libjst::insertion(store[4]), std::vector<alphabet_t>{});

    EXPECT_EQ(libjst::deletion(store[0]), 1u);
    EXPECT_EQ(libjst::deletion(store[1]), 1u);
    EXPECT_EQ(libjst::deletion(store[2]), this->insertion_sequence.size());
    EXPECT_EQ(libjst::deletion(store[3]), 0u);
    EXPECT_EQ(libjst::deletion(store[4]), 1u);

    EXPECT_RANGE_EQ(libjst::coverage(store[0]), (coverage_t{0, 0, 0, 1}));
    EXPECT_RANGE_EQ(libjst::coverage(store[1]), (coverage_t{1, 0, 0, 0}));
    EXPECT_RANGE_EQ(libjst::coverage(store[2]), (coverage_t{0, 0, 1, 0}));
    EXPECT_RANGE_EQ(libjst::coverage(store[3]), (coverage_t{0, 1, 0, 0}));
    EXPECT_RANGE_EQ(libjst::coverage(store[4]), (coverage_t{0, 0, 1, 1}));
}

TYPED_TEST(variant_store_covered_test, iterator)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    covered_store_t store{};

    EXPECT_EQ((store.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}}) - store.begin()), 0);
    EXPECT_EQ((store.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}}) - store.begin()), 2);
    EXPECT_EQ((store.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}}) - store.begin()), 1);
    EXPECT_EQ((store.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}}) - store.begin()), 4);

    auto it = store.begin();
    EXPECT_EQ(libjst::position(*it), 4u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_rank_to(3, alphabet_t{})}));
    EXPECT_EQ(libjst::deletion(*it), 1u);
    EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 0, 1}));

    ++it;
    EXPECT_EQ(libjst::position(*it), 112u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_rank_to(0, alphabet_t{})}));
    EXPECT_EQ(libjst::deletion(*it), 1u);
    EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{1, 0, 0, 0}));

    ++it;
    EXPECT_EQ(libjst::position(*it), 44u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), this->insertion_sequence);
    EXPECT_EQ(libjst::deletion(*it), this->insertion_sequence.size());
    EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 1, 0}));

    ++it;
    EXPECT_EQ(libjst::position(*it), 93u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), this->insertion_sequence);
    EXPECT_EQ(libjst::deletion(*it), 0u);
    EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 1, 0, 0}));

    ++it;
    EXPECT_EQ(libjst::position(*it), 154u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), std::vector<alphabet_t>{});
    EXPECT_EQ(libjst::deletion(*it), 1u);
    EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 1, 1}));

    ++it;
    EXPECT_EQ(it, store.end());
}

TYPED_TEST(variant_store_covered_test, serialise)
{
    using covered_store_t = typename TestFixture::covered_store_t;
    using coverage_t = typename TestFixture::coverage_t;
    using value_t = std::ranges::range_value_t<covered_store_t>;

    covered_store_t store_out{};
    covered_store_t store_in{};

    store_out.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}});
    store_out.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}});
    store_out.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}});
    store_out.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}});
    store_out.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}});

    std::stringstream archive_stream{};
    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        output_archive(store_out);
    }
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        input_archive(store_in);
    }

    for (size_t i = 0; i < std::ranges::size(store_out); ++i)
    {
        EXPECT_EQ(libjst::position(store_in[i]), libjst::position(store_out[i]));
        EXPECT_EQ(libjst::deletion(store_in[i]), libjst::deletion(store_out[i]));
        EXPECT_RANGE_EQ(libjst::insertion(store_in[i]), libjst::insertion(store_out[i]));
        EXPECT_RANGE_EQ(libjst::coverage(store_in[i]), libjst::coverage(store_out[i]));
    }
}
