// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <tuple>

// #include <cereal/archives/json.hpp>

#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/rcms/generic_delta.hpp>
#include <libjst/variant/breakpoint.hpp>

template <typename tuple_t>
struct generic_delta_test : ::testing::Test
{
    using source_type = std::tuple_element_t<0, tuple_t>;
    using coverage_type = std::tuple_element_t<1, tuple_t>;
    using breakpoint_type = libjst::breakpoint;
    using test_type = libjst::generic_delta<source_type, coverage_type>;

};

using test_types = ::testing::Types<std::tuple<std::string, std::vector<uint32_t>>
                                    >;
TYPED_TEST_SUITE(generic_delta_test, test_types);

TYPED_TEST(generic_delta_test, snv)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    test_type delta{breakpoint_type{9, 1}, source_type{"G"}, coverage_type{0, 2}};
    EXPECT_EQ(libjst::low_breakend(delta), 9);
    EXPECT_EQ(libjst::high_breakend(delta), 10);
    EXPECT_EQ(libjst::breakpoint_span(delta), 1);
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{"G"});
    EXPECT_RANGE_EQ(libjst::coverage(delta), (coverage_type{0, 2}));
}

TYPED_TEST(generic_delta_test, deletion)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    test_type delta{breakpoint_type{1, 7}, source_type{""}, coverage_type{1}};
    EXPECT_EQ(libjst::low_breakend(delta), 1);
    EXPECT_EQ(libjst::high_breakend(delta), 8);
    EXPECT_EQ(libjst::breakpoint_span(delta), 7);
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{""});
    EXPECT_RANGE_EQ(libjst::coverage(delta), coverage_type{1});
}

TYPED_TEST(generic_delta_test, insertion)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    test_type delta{breakpoint_type{13, 0}, source_type{"AAA"}, coverage_type{0, 1, 2, 3}};
    EXPECT_EQ(libjst::low_breakend(delta), 13);
    EXPECT_EQ(libjst::high_breakend(delta), 13);
    EXPECT_EQ(libjst::breakpoint_span(delta), 0);
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{"AAA"});
    EXPECT_RANGE_EQ(libjst::coverage(delta), (coverage_type{0, 1, 2, 3}));
}

TYPED_TEST(generic_delta_test, unbalanced_replacement)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    test_type delta{breakpoint_type{14, 3}, source_type{"A"}, coverage_type{0, 4}};
    EXPECT_EQ(libjst::low_breakend(delta), 14);
    EXPECT_EQ(libjst::high_breakend(delta), 17);
    EXPECT_EQ(libjst::breakpoint_span(delta), 3);
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{"A"});
    EXPECT_RANGE_EQ(libjst::coverage(delta), (coverage_type{0, 4}));
}

TYPED_TEST(generic_delta_test, assign)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    test_type delta{};
    EXPECT_EQ(libjst::low_breakend(delta), 0);
    EXPECT_EQ(libjst::high_breakend(delta), 0);
    EXPECT_EQ(libjst::breakpoint_span(delta), 0);
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{""});
    EXPECT_RANGE_EQ(libjst::coverage(delta), (coverage_type{}));

    libjst::get_breakpoint(delta) = breakpoint_type{2, 4};
    EXPECT_EQ(libjst::low_breakend(delta), 2);
    EXPECT_EQ(libjst::high_breakend(delta), 6);

    libjst::alt_sequence(delta) = source_type{"AAA"};
    EXPECT_RANGE_EQ(libjst::alt_sequence(delta), source_type{"AAA"});

    libjst::coverage(delta) = coverage_type{0, 1, 2};
    EXPECT_RANGE_EQ(libjst::coverage(delta), (coverage_type{0, 1, 2}));
}

TYPED_TEST(generic_delta_test, lvalue_reference)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    EXPECT_TRUE((std::same_as<decltype(libjst::get_breakpoint(std::declval<test_type &>())), breakpoint_type &>));
    EXPECT_TRUE((std::same_as<decltype(libjst::alt_sequence(std::declval<test_type &>())), source_type &>));
    EXPECT_TRUE((std::same_as<decltype(libjst::coverage(std::declval<test_type &>())), coverage_type &>));
}

TYPED_TEST(generic_delta_test, lvalue_const_reference)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    EXPECT_TRUE((std::same_as<decltype(libjst::get_breakpoint(std::declval<test_type const &>())), breakpoint_type const &>));
    EXPECT_TRUE((std::same_as<decltype(libjst::alt_sequence(std::declval<test_type const &>())), source_type const &>));
    EXPECT_TRUE((std::same_as<decltype(libjst::coverage(std::declval<test_type const &>())), coverage_type const &>));
}

TYPED_TEST(generic_delta_test, rvalue_reference)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    EXPECT_TRUE((std::same_as<decltype(libjst::get_breakpoint(std::declval<test_type &&>())), breakpoint_type &&>));
    EXPECT_TRUE((std::same_as<decltype(libjst::alt_sequence(std::declval<test_type &&>())), source_type &&>));
    EXPECT_TRUE((std::same_as<decltype(libjst::coverage(std::declval<test_type &&>())), coverage_type &&>));
}

TYPED_TEST(generic_delta_test, const_rvalue_reference)
{
    using test_type = typename TestFixture::test_type;
    using coverage_type = typename TestFixture::coverage_type;
    using source_type = typename TestFixture::source_type;
    using breakpoint_type = typename TestFixture::breakpoint_type;

    EXPECT_TRUE((std::same_as<decltype(libjst::get_breakpoint(std::declval<test_type const &&>())), breakpoint_type const &&>));
    EXPECT_TRUE((std::same_as<decltype(libjst::alt_sequence(std::declval<test_type const &&>())), source_type const &&>));
    EXPECT_TRUE((std::same_as<decltype(libjst::coverage(std::declval<test_type const &&>())), coverage_type const &&>));
}


// TYPED_TEST(generic_delta_test, serialise)
// {
//     using generic_variant_t = typename TestFixture::generic_variant_t;
//     using alphabet_t = typename TestFixture::alphabet_t;

//     generic_variant_t var_sub_out{0, this->insertion_sequence, (uint32_t)std::ranges::size(this->insertion_sequence)};
//     generic_variant_t var_del_out{1234, std::vector<alphabet_t>{}, 15};
//     generic_variant_t var_ins_out{((1 << 30) - 1), this->insertion_sequence, 0};

//     generic_variant_t var_sub_in{};
//     generic_variant_t var_ins_in{};
//     generic_variant_t var_del_in{};

//     std::stringstream archive_stream{};
//     {
//         cereal::JSONOutputArchive output_archive(archive_stream);
//         output_archive(var_sub_out);
//         output_archive(var_del_out);
//         output_archive(var_ins_out);
//     }
//     {
//         cereal::JSONInputArchive input_archive(archive_stream);
//         input_archive(var_sub_in);
//         input_archive(var_del_in);
//         input_archive(var_ins_in);
//     }

//     EXPECT_EQ(libjst::position(var_sub_in), libjst::position(var_sub_out));
//     EXPECT_EQ(libjst::deletion(var_sub_in), libjst::deletion(var_sub_out));
//     EXPECT_RANGE_EQ(libjst::insertion(var_sub_in), libjst::insertion(var_sub_out));

//     EXPECT_EQ(libjst::position(var_del_in), libjst::position(var_del_out));
//     EXPECT_EQ(libjst::deletion(var_del_in), libjst::deletion(var_del_out));
//     EXPECT_RANGE_EQ(libjst::insertion(var_del_in), libjst::insertion(var_del_out));

//     EXPECT_EQ(libjst::position(var_ins_in), libjst::position(var_ins_out));
//     EXPECT_EQ(libjst::deletion(var_ins_in), libjst::deletion(var_ins_out));
//     EXPECT_RANGE_EQ(libjst::insertion(var_ins_in), libjst::insertion(var_ins_out));
// }
