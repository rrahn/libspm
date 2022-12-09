// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>

// #include <cereal/types/vector.hpp>
// #include <cereal/archives/json.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

#include <libcontrib/seqan/alphabet.hpp>

#include <libjst/variant/single_base_replacement_store.hpp>

template <typename alphabet_t>
struct single_base_replacement_test : public ::testing::Test
{
    using alphabet_type = alphabet_t;
    using store_type = libjst::single_base_replacement_store<alphabet_type>;

    alphabet_type value0{seqan3::assign_char_to('A', alphabet_t{})};
    alphabet_type value1{seqan3::assign_char_to('T', alphabet_t{})};

    auto to_sequence(char const c) {
        return std::array{seqan3::assign_char_to(c, alphabet_t{})};
    }
};

using test_types = ::testing::Types<jst::contrib::dna4,
                                    seqan3::dna4,
                                    jst::contrib::dna5,
                                    jst::contrib::dna15
                                >;
TYPED_TEST_SUITE(single_base_replacement_test, test_types);

TYPED_TEST(single_base_replacement_test, construction)
{
    using store_t = typename TestFixture::store_type;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<store_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<store_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<store_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<store_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<store_t>);
    EXPECT_TRUE(std::is_destructible_v<store_t>);
}

TYPED_TEST(single_base_replacement_test, range_concept)
{
    using store_t = typename TestFixture::store_type;
    EXPECT_TRUE(std::ranges::random_access_range<store_t>);
    EXPECT_TRUE(std::ranges::sized_range<store_t>);
}

TYPED_TEST(single_base_replacement_test, value_concept)
{
    using store_t = typename TestFixture::store_type;
    using value_t = std::ranges::range_value_t<store_t>;

    // EXPECT_TRUE(libjst::sequence_alternative<value_t>);
    EXPECT_EQ(libjst::breakpoint_span(value_t{this->value0}), 1ul);
    EXPECT_RANGE_EQ(libjst::alt_sequence(value_t{this->value0}), this->to_sequence('A'));
    EXPECT_EQ(libjst::effective_size(value_t{this->value0}), 0);
    EXPECT_EQ(libjst::alt_kind(value_t{this->value0}), libjst::alternate_sequence_kind::replacement);
}

TYPED_TEST(single_base_replacement_test, reserve)
{
    using store_t = typename TestFixture::store_type;

    store_t store{};

    std::size_t old_cap = store.capacity();
    store.reserve(old_cap + 1);
    EXPECT_LT(old_cap, store.capacity());
}

TYPED_TEST(single_base_replacement_test, resize)
{
    using store_t = typename TestFixture::store_type;

    store_t store{};

    std::size_t old_size = store.size();
    store.resize(old_size + 1);
    EXPECT_EQ(old_size + 1, store.size());
}

TYPED_TEST(single_base_replacement_test, push_back)
{
    using store_t = typename TestFixture::store_type;

    store_t store{};

    EXPECT_TRUE(std::ranges::empty(store));
    store.push_back(this->value0);
    store.push_back(this->value1);

    EXPECT_FALSE(std::ranges::empty(store));
    EXPECT_RANGE_EQ(libjst::alt_sequence(store[0]), this->to_sequence('A'));
    EXPECT_RANGE_EQ(libjst::alt_sequence(store[1]), this->to_sequence('T'));
}

// TYPED_TEST(single_base_replacement_test, emplace)
// {
//     using store_type = typename TestFixture::store_type;
//     using coverage_t = typename TestFixture::coverage_t;

//     store_type store{};

//     EXPECT_EQ((store.emplace(this->snp0, coverage_t{0, 0, 0, 1}) - store.begin()), 0);
//     EXPECT_EQ((store.emplace(this->var0, coverage_t{0, 0, 1, 0}) - store.begin()), 1);
//     EXPECT_EQ((store.emplace(this->var1, coverage_t{0, 1, 0, 0}) - store.begin()), 2);
//     EXPECT_EQ((store.emplace(this->snp1, coverage_t{1, 0, 0, 0}) - store.begin()), 1);
//     EXPECT_EQ((store.emplace(this->var2, coverage_t{0, 0, 1, 1}) - store.begin()), 4);
// }

// TYPED_TEST(single_base_replacement_test, size)
// {
//     using store_type = typename TestFixture::store_type;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<store_type>;

//     store_type store{};
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

// TYPED_TEST(single_base_replacement_test, subscript)
// {
//     using store_type = typename TestFixture::store_type;
//     using alphabet_t = typename TestFixture::alphabet_t;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<store_type>;

//     store_type store{};

//     EXPECT_EQ((store.insert(value_t{this->snp0, coverage_t{0, 0, 0, 1}}) - store.begin()), 0);
//     EXPECT_EQ((store.insert(value_t{this->var0, coverage_t{0, 0, 1, 0}}) - store.begin()), 1);
//     EXPECT_EQ((store.insert(value_t{this->var1, coverage_t{0, 1, 0, 0}}) - store.begin()), 2);
//     EXPECT_EQ((store.insert(value_t{this->snp1, coverage_t{1, 0, 0, 0}}) - store.begin()), 1);
//     EXPECT_EQ((store.insert(value_t{this->var2, coverage_t{0, 0, 1, 1}}) - store.begin()), 4);

//     EXPECT_EQ(libjst::position(store[0]), 4u);
//     EXPECT_EQ(libjst::position(store[1]), 112u);
//     EXPECT_EQ(libjst::position(store[2]), 44u);
//     EXPECT_EQ(libjst::position(store[3]), 93u);
//     EXPECT_EQ(libjst::position(store[4]), 154u);

//     EXPECT_RANGE_EQ(libjst::insertion(store[0]), (std::vector{seqan3::assign_char_to('T', alphabet_t{})}));
//     EXPECT_RANGE_EQ(libjst::insertion(store[1]), (std::vector{seqan3::assign_char_to('A', alphabet_t{})}));
//     EXPECT_RANGE_EQ(libjst::insertion(store[2]), this->insertion_sequence);
//     EXPECT_RANGE_EQ(libjst::insertion(store[3]), this->insertion_sequence);
//     EXPECT_RANGE_EQ(libjst::insertion(store[4]), std::vector<alphabet_t>{});

//     EXPECT_EQ(libjst::deletion(store[0]), 1u);
//     EXPECT_EQ(libjst::deletion(store[1]), 1u);
//     EXPECT_EQ(libjst::deletion(store[2]), this->insertion_sequence.size());
//     EXPECT_EQ(libjst::deletion(store[3]), 0u);
//     EXPECT_EQ(libjst::deletion(store[4]), 1u);

//     EXPECT_RANGE_EQ(libjst::coverage(store[0]), (coverage_t{0, 0, 0, 1}));
//     EXPECT_RANGE_EQ(libjst::coverage(store[1]), (coverage_t{1, 0, 0, 0}));
//     EXPECT_RANGE_EQ(libjst::coverage(store[2]), (coverage_t{0, 0, 1, 0}));
//     EXPECT_RANGE_EQ(libjst::coverage(store[3]), (coverage_t{0, 1, 0, 0}));
//     EXPECT_RANGE_EQ(libjst::coverage(store[4]), (coverage_t{0, 0, 1, 1}));
// }

// TYPED_TEST(single_base_replacement_test, iterator)
// {
//     using store_type = typename TestFixture::store_type;
//     using alphabet_t = typename TestFixture::alphabet_t;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<store_type>;

//     store_type store{};

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

// TYPED_TEST(single_base_replacement_test, serialise)
// {
//     using store_type = typename TestFixture::store_type;
//     using coverage_t = typename TestFixture::coverage_t;
//     using value_t = std::ranges::range_value_t<store_type>;

//     store_type store_out{};
//     store_type store_in{};

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
