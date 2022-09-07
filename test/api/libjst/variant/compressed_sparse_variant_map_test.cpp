// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// #include <cereal/types/vector.hpp>
// #include <cereal/archives/json.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/utility/bit_vector.hpp>
#include <libjst/variant/single_base_replacement_store.hpp>
#include <libjst/variant/compressed_sparse_variant_map.hpp>

template <typename alphabet_t>
struct compressed_sparse_variant_map_test : public ::testing::Test
{
    using alphabet_type = alphabet_t;
    using alt_store_type = libjst::single_base_replacement_store<alphabet_type>;
    using alt_value_type = std::ranges::range_value_t<alt_store_type>;
    using coverage_t = libjst::bit_vector<>;

    using map_type = libjst::compressed_sparse_variant_map<alt_store_type, coverage_t>;
    using test_variant_type = std::tuple<libjst::breakpoint, alt_value_type, coverage_t>;

    test_variant_type test_snvA{libjst::breakpoint{10}, alt_value_type{seqan3::assign_char_to('A', alphabet_t{})}, coverage_t{1, 0, 0, 0}};
    test_variant_type test_snvC{libjst::breakpoint{15}, alt_value_type{seqan3::assign_char_to('C', alphabet_t{})}, coverage_t{0, 1, 0, 0}};
    test_variant_type test_snvG{libjst::breakpoint{10}, alt_value_type{seqan3::assign_char_to('G', alphabet_t{})}, coverage_t{0, 0, 1, 0}};
    test_variant_type test_snvT{libjst::breakpoint{7},  alt_value_type{seqan3::assign_char_to('T', alphabet_t{})}, coverage_t{0, 0, 0, 1}};

    auto to_sequence(char const c) {
        return std::array{seqan3::assign_char_to(c, alphabet_t{})};
    }
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4,
                                    jst::contrib::dna5,
                                    jst::contrib::dna15
                                >;
TYPED_TEST_SUITE(compressed_sparse_variant_map_test, test_types);

TYPED_TEST(compressed_sparse_variant_map_test, construction)
{
    using map_t = typename TestFixture::map_type;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<map_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<map_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<map_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<map_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<map_t>);
    EXPECT_TRUE(std::is_destructible_v<map_t>);
}

TYPED_TEST(compressed_sparse_variant_map_test, range_concept)
{
    using map_t = typename TestFixture::map_type;
    EXPECT_TRUE(std::ranges::random_access_range<map_t>);
    EXPECT_TRUE(std::ranges::sized_range<map_t>);
}

TYPED_TEST(compressed_sparse_variant_map_test, iterator_concept)
{
    using map_t = typename TestFixture::map_type;
    using iterator = std::ranges::iterator_t<map_t>;
    using const_iterator = std::ranges::iterator_t<map_t const>;

    EXPECT_TRUE((std::constructible_from<const_iterator, iterator>));
}

TYPED_TEST(compressed_sparse_variant_map_test, value_concept)
{
    // using map_t = typename TestFixture::map_type;
    // using value_t = std::ranges::range_value_t<map_t>;

    // EXPECT_TRUE(libjst::sequence_alternative<value_t>);
    // EXPECT_TRUE(libjst::sequence_variant<value_t>);
    // EXPECT_TRUE(libjst::covered_sequence_variant<value_t>);
}

TYPED_TEST(compressed_sparse_variant_map_test, insert)
{
    using map_t = typename TestFixture::map_type;

    map_t map{};

    EXPECT_TRUE(std::ranges::empty(map));

    auto map_insert = [&] (auto &&...args) {
        return map.emplace(std::forward<decltype(args)>(args)...);
    };

    auto pos_snvC = std::distance(map.begin(), std::apply(map_insert, this->test_snvC));
    EXPECT_EQ(pos_snvC, 0);
    auto pos_snvG = std::distance(map.begin(), std::apply(map_insert, this->test_snvG));
    EXPECT_EQ(pos_snvG, 0);
    auto pos_snvA = std::distance(map.begin(), std::apply(map_insert, this->test_snvA));
    EXPECT_EQ(pos_snvA, 0);
    auto pos_snvT = std::distance(map.begin(), std::apply(map_insert, this->test_snvT));
    EXPECT_EQ(pos_snvT, 0);

    auto it = map.begin();
    EXPECT_EQ(libjst::position(*it), get<0>(this->test_snvT));
    EXPECT_RANGE_EQ(libjst::alt_sequence(*it), this->to_sequence('T'));
    ++it;
    EXPECT_EQ(libjst::position(*it), get<0>(this->test_snvA));
    EXPECT_RANGE_EQ(libjst::alt_sequence(*it), this->to_sequence('A'));
    ++it;
    EXPECT_EQ(libjst::position(*it), get<0>(this->test_snvG));
    EXPECT_RANGE_EQ(libjst::alt_sequence(*it), this->to_sequence('G'));
    ++it;
    EXPECT_EQ(libjst::position(*it), get<0>(this->test_snvC));
    EXPECT_RANGE_EQ(libjst::alt_sequence(*it), this->to_sequence('C'));
}

// TYPED_TEST(compressed_sparse_variant_map_test, emplace)
// {
//     using map_type = typename TestFixture::map_type;
//     using coverage_t = typename TestFixture::coverage_t;

//     map_type store{};

//     EXPECT_EQ((store.emplace(this->snp0, coverage_t{0, 0, 0, 1}) - store.begin()), 0);
//     EXPECT_EQ((store.emplace(this->var0, coverage_t{0, 0, 1, 0}) - store.begin()), 1);
//     EXPECT_EQ((store.emplace(this->var1, coverage_t{0, 1, 0, 0}) - store.begin()), 2);
//     EXPECT_EQ((store.emplace(this->snp1, coverage_t{1, 0, 0, 0}) - store.begin()), 1);
//     EXPECT_EQ((store.emplace(this->var2, coverage_t{0, 0, 1, 1}) - store.begin()), 4);
// }

// TYPED_TEST(compressed_sparse_variant_map_test, size)
// {
//     using map_type = typename TestFixture::map_type;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<map_type>;

//     map_type store{};
//     coverage_t cov{0, 1, 0, 1};
//     EXPECT_EQ(store.size(), 0u);
//     EXPECT_NO_THROW((store.insert(value_t{this->snp0, cov})));
//     EXPECT_EQ(store.size(), 1u);
//     EXPECT_NO_THROW((store.insert(value_t{this->var0, cov})));
//     EXPECT_EQ(store.size(), 2u);
//     EXPECT_NO_THROW((store.insert(value_t{this->var1, cov})));
//     EXPECT_EQ(store.size(), 3u);
//     EXPECT_NO_THROW((store.insert(value_t{this->snp1, cov})));
//     EXPECT_EQ(store.size(), 4u);
//     EXPECT_NO_THROW((store.insert(value_t{this->var2, cov})));
//     EXPECT_EQ(store.size(), 5u);
// }

// TYPED_TEST(compressed_sparse_variant_map_test, iterator)
// {
//     using map_type = typename TestFixture::map_type;
//     using alphabet_t = typename TestFixture::alphabet_t;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<map_type>;

//     map_type store{};

//     EXPECT_EQ((store.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}}) - store.begin()), 0);
//     EXPECT_EQ((store.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}}) - store.begin()), 1);
//     EXPECT_EQ((store.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}}) - store.begin()), 2);
//     EXPECT_EQ((store.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}}) - store.begin()), 1);
//     EXPECT_EQ((store.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}}) - store.begin()), 4);

//     auto it = store.begin();
//     EXPECT_EQ(libjst::position(*it), 4u);
//     EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_char_to('T', alphabet_t{})}));
//     EXPECT_EQ(libjst::deletion(*it), 1u);
//     EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 0, 1}));

//     ++it;
//     EXPECT_EQ(libjst::position(*it), 112u);
//     EXPECT_RANGE_EQ(libjst::insertion(*it), (std::vector{seqan3::assign_char_to('A', alphabet_t{})}));
//     EXPECT_EQ(libjst::deletion(*it), 1u);
//     EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{1, 0, 0, 0}));

//     ++it;
//     EXPECT_EQ(libjst::position(*it), 44u);
//     EXPECT_RANGE_EQ(libjst::insertion(*it), this->insertion_sequence);
//     EXPECT_EQ(libjst::deletion(*it), this->insertion_sequence.size());
//     EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 1, 0}));

//     ++it;
//     EXPECT_EQ(libjst::position(*it), 93u);
//     EXPECT_RANGE_EQ(libjst::insertion(*it), this->insertion_sequence);
//     EXPECT_EQ(libjst::deletion(*it), 0u);
//     EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 1, 0, 0}));

//     ++it;
//     EXPECT_EQ(libjst::position(*it), 154u);
//     EXPECT_RANGE_EQ(libjst::insertion(*it), std::vector<alphabet_t>{});
//     EXPECT_EQ(libjst::deletion(*it), 1u);
//     EXPECT_RANGE_EQ(libjst::coverage(*it), (coverage_t{0, 0, 1, 1}));

//     ++it;
//     EXPECT_EQ(it, store.end());
// }

// TYPED_TEST(compressed_sparse_variant_map_test, serialise)
// {
//     using map_type = typename TestFixture::map_type;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<map_type>;

//     map_type store_out{};
//     map_type store_in{};

//     store_out.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}});
//     store_out.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}});
//     store_out.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}});
//     store_out.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}});
//     store_out.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}});

//     std::stringstream archive_stream{};
//     {
//         cereal::JSONOutputArchive output_archive(archive_stream);
//         output_archive(store_out);
//     }
//     {
//         cereal::JSONInputArchive input_archive(archive_stream);
//         input_archive(store_in);
//     }

//     for (size_t i = 0; i < std::ranges::size(store_out); ++i)
//     {
//         EXPECT_EQ(libjst::position(store_in[i]), libjst::position(store_out[i]));
//         EXPECT_EQ(libjst::deletion(store_in[i]), libjst::deletion(store_out[i]));
//         EXPECT_RANGE_EQ(libjst::insertion(store_in[i]), libjst::insertion(store_out[i]));
//         EXPECT_RANGE_EQ(libjst::coverage(store_in[i]), libjst::coverage(store_out[i]));
//     }
// }
