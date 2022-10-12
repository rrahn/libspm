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

#include <libjst/sequence_variant/variant_snp.hpp>
#include <libjst/sequence_variant/variant_generic.hpp>
#include <libjst/sequence_variant/variant_store_composite.hpp>

template <typename alphabet_type>
struct variant_store_composite_test : public ::testing::Test
{
    using alphabet_t = alphabet_type;
    using snp_variant_t = libjst::snp_variant<alphabet_t>;
    using generic_variant_t = libjst::generic_variant<alphabet_t>;

    using snp_store_t = std::vector<snp_variant_t>;
    using generic_store_t = std::vector<generic_variant_t>;
    using composite_store_t = libjst::variant_store_composite<snp_store_t, generic_store_t>;

    auto make_insertion()
    {
        return seqan3::test::generate_sequence<alphabet_t>(10);
    }
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4,
                                    jst::contrib::dna5,
                                    jst::contrib::dna15
                                >;
TYPED_TEST_SUITE(variant_store_composite_test, test_types);

TYPED_TEST(variant_store_composite_test, construction)
{
    using composite_store_t = typename TestFixture::composite_store_t;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<composite_store_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<composite_store_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<composite_store_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<composite_store_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<composite_store_t>);
    EXPECT_TRUE(std::is_destructible_v<composite_store_t>);
}

TYPED_TEST(variant_store_composite_test, concept)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    EXPECT_TRUE(std::ranges::random_access_range<composite_store_t>);
    EXPECT_TRUE(libjst::sequence_variant_store<composite_store_t>);
    EXPECT_TRUE(libjst::sequence_variant_store<composite_store_t &>);
    EXPECT_TRUE(libjst::sequence_variant_store<composite_store_t const &>);
    EXPECT_FALSE(libjst::covered_sequence_variant_store<composite_store_t>);
    EXPECT_FALSE(libjst::covered_sequence_variant_store<composite_store_t &>);
    EXPECT_FALSE(libjst::covered_sequence_variant_store<composite_store_t const &>);
}

TYPED_TEST(variant_store_composite_test, type_traits)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;
    using value_t = std::ranges::range_value_t<composite_store_t>;
    using reference_t = std::ranges::range_reference_t<composite_store_t const>;

    // constructibility of the value type
    EXPECT_TRUE((std::constructible_from<value_t, value_t &>));
    EXPECT_TRUE((std::constructible_from<value_t, value_t const &>));
    EXPECT_TRUE((std::constructible_from<value_t, value_t>));
    EXPECT_TRUE((std::constructible_from<value_t, value_t &&>));
    EXPECT_TRUE((std::constructible_from<value_t, value_t const &&>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t &>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t const &>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t &&>));
    EXPECT_TRUE((std::constructible_from<value_t, reference_t const &&>));
    EXPECT_TRUE((std::constructible_from<value_t, snp_variant_t &>));
    EXPECT_TRUE((std::constructible_from<value_t, snp_variant_t const &>));
    EXPECT_TRUE((std::constructible_from<value_t, snp_variant_t>));
    EXPECT_TRUE((std::constructible_from<value_t, snp_variant_t &&>));
    EXPECT_TRUE((std::constructible_from<value_t, snp_variant_t const &&>));
    EXPECT_TRUE((std::constructible_from<value_t, generic_variant_t &>));
    EXPECT_TRUE((std::constructible_from<value_t, generic_variant_t const &>));
    EXPECT_TRUE((std::constructible_from<value_t, generic_variant_t>));
    EXPECT_TRUE((std::constructible_from<value_t, generic_variant_t &&>));
    EXPECT_TRUE((std::constructible_from<value_t, generic_variant_t const &&>));

    // constructibility of the reference type
    EXPECT_TRUE((std::constructible_from<reference_t, reference_t &>));
    EXPECT_TRUE((std::constructible_from<reference_t, reference_t const &>));
    EXPECT_TRUE((std::constructible_from<reference_t, reference_t>));
    EXPECT_TRUE((std::constructible_from<reference_t, reference_t &&>));
    EXPECT_TRUE((std::constructible_from<reference_t, reference_t const &&>));
    EXPECT_TRUE((std::constructible_from<reference_t, value_t &>));
    EXPECT_TRUE((std::constructible_from<reference_t, value_t const &>));
    // Constructible from rvalue value type since we store const & types inside.
    EXPECT_TRUE((std::constructible_from<reference_t, value_t>));
    EXPECT_TRUE((std::constructible_from<reference_t, value_t &&>));
    EXPECT_TRUE((std::constructible_from<reference_t, value_t const &&>));
    EXPECT_TRUE((std::constructible_from<reference_t, snp_variant_t &>));
    EXPECT_TRUE((std::constructible_from<reference_t, snp_variant_t const &>));
    EXPECT_FALSE((std::constructible_from<reference_t, snp_variant_t>));
    EXPECT_FALSE((std::constructible_from<reference_t, snp_variant_t &&>));
    EXPECT_FALSE((std::constructible_from<reference_t, snp_variant_t const &&>));
    EXPECT_TRUE((std::constructible_from<reference_t, generic_variant_t &>));
    EXPECT_TRUE((std::constructible_from<reference_t, generic_variant_t const &>));
    EXPECT_FALSE((std::constructible_from<reference_t, generic_variant_t>));
    EXPECT_FALSE((std::constructible_from<reference_t, generic_variant_t &&>));
    EXPECT_FALSE((std::constructible_from<reference_t, generic_variant_t const &&>));

    EXPECT_EQ(sizeof(reference_t), 16u);

    snp_variant_t snp{4, seqan3::assign_char_to('T', alphabet_t{})};
    reference_t any_ref{snp};
    value_t any_val{snp_variant_t{snp}};
    EXPECT_EQ(libjst::position(snp), 4u);
    EXPECT_EQ(libjst::position(any_ref), 4u);
    EXPECT_EQ(libjst::position(any_val), 4u);

    snp = snp_variant_t{10, seqan3::assign_char_to('A', alphabet_t{})};
    EXPECT_EQ(libjst::position(snp), 10u);
    EXPECT_EQ(libjst::position(any_ref), 10u);
    EXPECT_EQ(libjst::position(any_val), 4u);
}

TYPED_TEST(variant_store_composite_test, insert)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store{};
    EXPECT_EQ((store.insert(snp_variant_t{4, seqan3::assign_char_to('T', alphabet_t{})}) - store.begin()), 0);
    EXPECT_EQ((store.insert(generic_variant_t{44, this->make_insertion(), 10}) - store.begin()), 1);
    EXPECT_EQ((store.insert(generic_variant_t{93, this->make_insertion(), 0}) - store.begin()), 2);
    EXPECT_EQ((store.insert(snp_variant_t{112, seqan3::assign_char_to('A', alphabet_t{})}) - store.begin()), 1);
    EXPECT_EQ((store.insert(generic_variant_t{154, {}, 1}) - store.begin()), 4);
}

TYPED_TEST(variant_store_composite_test, emplace)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    // using snp_variant_t = typename TestFixture::snp_variant_t;
    // using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store{};
    EXPECT_EQ((store.emplace(uint32_t{4}, seqan3::assign_char_to('T', alphabet_t{})) - store.begin()), 0);
    EXPECT_EQ((store.emplace(uint32_t{44}, this->make_insertion(), uint32_t{10}) - store.begin()), 1);
    EXPECT_EQ((store.emplace(uint32_t{93}, this->make_insertion(), uint32_t{0}) - store.begin()), 2);
    EXPECT_EQ((store.emplace(uint32_t{112}, seqan3::assign_char_to('A', alphabet_t{})) - store.begin()), 1);
    EXPECT_EQ((store.emplace(uint32_t{154}, std::vector<alphabet_t>{}, uint32_t{1}) - store.begin()), 4);
}

TYPED_TEST(variant_store_composite_test, size)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store{};
    EXPECT_EQ(store.size(), 0u);
    EXPECT_NO_THROW((store.insert(snp_variant_t{4, seqan3::assign_char_to('T', alphabet_t{})})));
    EXPECT_EQ(store.size(), 1u);
    EXPECT_NO_THROW((store.insert(generic_variant_t{44, this->make_insertion(), 10})));
    EXPECT_EQ(store.size(), 2u);
    EXPECT_NO_THROW((store.insert(generic_variant_t{93, this->make_insertion(), 0})));
    EXPECT_EQ(store.size(), 3u);
    EXPECT_NO_THROW((store.insert(snp_variant_t{112, seqan3::assign_char_to('A', alphabet_t{})})));
    EXPECT_EQ(store.size(), 4u);
    EXPECT_NO_THROW((store.insert(generic_variant_t{154, {}, 1})));
    EXPECT_EQ(store.size(), 5u);
}

TYPED_TEST(variant_store_composite_test, subscript)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store{};
    auto ins = this->make_insertion();

    EXPECT_NO_THROW((store.insert(snp_variant_t{4, seqan3::assign_char_to('T', alphabet_t{})})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{44, ins, (uint32_t) ins.size()})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{93, ins, 0})));
    EXPECT_NO_THROW((store.insert(snp_variant_t{112, seqan3::assign_char_to('A', alphabet_t{})})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{154, {}, 1})));

    EXPECT_EQ(libjst::position(store[0]), 4u);
    EXPECT_EQ(libjst::position(store[1]), 112u);
    EXPECT_EQ(libjst::position(store[2]), 44u);
    EXPECT_EQ(libjst::position(store[3]), 93u);
    EXPECT_EQ(libjst::position(store[4]), 154u);

    EXPECT_RANGE_EQ(libjst::insertion(store[0]), (std::vector{seqan3::assign_char_to('T', alphabet_t{})}));
    EXPECT_RANGE_EQ(libjst::insertion(store[1]), (std::vector{seqan3::assign_char_to('A', alphabet_t{})}));
    EXPECT_RANGE_EQ(libjst::insertion(store[2]), ins);
    EXPECT_RANGE_EQ(libjst::insertion(store[3]), ins);
    EXPECT_RANGE_EQ(libjst::insertion(store[4]), std::vector<alphabet_t>{});

    EXPECT_EQ(libjst::deletion(store[0]), 1u);
    EXPECT_EQ(libjst::deletion(store[1]), 1u);
    EXPECT_EQ(libjst::deletion(store[2]), ins.size());
    EXPECT_EQ(libjst::deletion(store[3]), 0u);
    EXPECT_EQ(libjst::deletion(store[4]), 1u);
}

TYPED_TEST(variant_store_composite_test, iterator)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store{};
    auto ins = this->make_insertion();

    EXPECT_NO_THROW((store.insert(snp_variant_t{4, seqan3::assign_char_to('T', alphabet_t{})})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{44, ins, (uint32_t) ins.size()})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{93, ins, 0})));
    EXPECT_NO_THROW((store.insert(snp_variant_t{112, seqan3::assign_char_to('A', alphabet_t{})})));
    EXPECT_NO_THROW((store.insert(generic_variant_t{154, {}, 1})));

    auto it = store.begin();
    EXPECT_EQ(libjst::position(*it), 4u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_char_to('T', alphabet_t{})}));
    EXPECT_EQ(libjst::deletion(*it), 1u);

    ++it;
    EXPECT_EQ(libjst::position(*it), 112u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_char_to('A', alphabet_t{})}));
    EXPECT_EQ(libjst::deletion(*it), 1u);

    ++it;
    EXPECT_EQ(libjst::position(*it), 44u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), ins);
    EXPECT_EQ(libjst::deletion(*it), ins.size());

    ++it;
    EXPECT_EQ(libjst::position(*it), 93u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), ins);
    EXPECT_EQ(libjst::deletion(*it), 0u);

    ++it;
    EXPECT_EQ(libjst::position(*it), 154u);
    EXPECT_RANGE_EQ(libjst::insertion(*it), std::vector<alphabet_t>{});
    EXPECT_EQ(libjst::deletion(*it), 1u);

    ++it;
    EXPECT_EQ(it, store.end());
}

TYPED_TEST(variant_store_composite_test, serialise)
{
    using composite_store_t = typename TestFixture::composite_store_t;
    using alphabet_t = typename TestFixture::alphabet_t;
    using snp_variant_t = typename TestFixture::snp_variant_t;
    using generic_variant_t = typename TestFixture::generic_variant_t;

    composite_store_t store_out{};
    composite_store_t store_in{};
    auto ins = this->make_insertion();

    EXPECT_NO_THROW((store_out.insert(snp_variant_t{4, seqan3::assign_char_to('T', alphabet_t{})})));
    EXPECT_NO_THROW((store_out.insert(generic_variant_t{44, ins, (uint32_t) ins.size()})));
    EXPECT_NO_THROW((store_out.insert(generic_variant_t{93, ins, 0})));
    EXPECT_NO_THROW((store_out.insert(snp_variant_t{112, seqan3::assign_char_to('A', alphabet_t{})})));
    EXPECT_NO_THROW((store_out.insert(generic_variant_t{154, {}, 1})));

    generic_variant_t var_sub_in{};
    generic_variant_t var_ins_in{};
    generic_variant_t var_del_in{};

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
    }
}
