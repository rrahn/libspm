// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>
#include <string>
#include <type_traits>

#include <seqan3/alphabet/adaptation/char.hpp>

#include <libcontrib/type_traits.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/container/journaled_sequence_group.hpp>
#include <libjst/variant/variant_generic.hpp>
#include <libjst/utility/bit_vector.hpp>

template <typename variant_t>
class covered_variant
{
private:
    using coverage_t = libjst::bit_vector<>;

    variant_t _variant{};
    coverage_t _coverage{};

public:

    using variant_type = variant_t;
    using coverage_type = coverage_t;

    covered_variant(variant_t variant, coverage_t coverage) :
        _variant{std::move(variant)},
        _coverage{std::move(coverage)}
    {}

private:

    template <typename this_t, typename member_t>
    using fwd_t = std::conditional_t<std::is_lvalue_reference_v<member_t>,
                                        member_t,
                                        jst::contrib::member_type_t<this_t, std::remove_cvref_t<member_t>>>;

    template <typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, covered_variant>
    constexpr friend auto tag_invoke(std::tag_t<libjst::coverage>, this_t &&me) noexcept
        -> fwd_t<this_t, coverage_t>
    {
        return (fwd_t<this_t, coverage_t> &&)me._coverage;
    }

    template <typename cpo_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, covered_variant> &&
                    std::tag_invocable<cpo_t, fwd_t<this_t, variant_t>>
    constexpr friend auto tag_invoke(cpo_t cpo, this_t &&me)
        noexcept(std::is_nothrow_tag_invocable_v<cpo_t, fwd_t<this_t, variant_t>>)
        -> std::tag_invoke_result_t<cpo_t, fwd_t<this_t, variant_t>>
    {
        return std::tag_invoke(cpo, (fwd_t<this_t, variant_t> &&)me._variant);
    }
};

template <typename t>
struct journaled_sequence_group_fixture : public ::testing::Test
{
    using source_t = typename t::first_type;
    using store_t = typename t::second_type;

    using type = libjst::journaled_sequence_group<std::views::all_t<source_t const &>, store_t>;

    source_t _test_source{"test source sequence"};

    source_t const & source() const noexcept
    {
        return _test_source;
    }

    store_t generate_variants() const noexcept
    {
        using value_t = std::ranges::range_value_t<store_t>;

        using variant_t = typename value_t::variant_type;
        using coverage_t = typename value_t::coverage_type;

        store_t tmp{};
        tmp.emplace_back(variant_t{4, std::vector<char>{'S', 'U', 'B'}, 3}, coverage_t{0, 0, 0, 1});
        tmp.emplace_back(variant_t{9, std::vector<char>{}, 2}, coverage_t{1, 1, 0, 0});
        return tmp;
    }
};

using test_types = ::testing::Types<
    std::pair<std::string, std::vector<covered_variant<libjst::generic_variant<char>>>>
>;
TYPED_TEST_SUITE(journaled_sequence_group_fixture, test_types, );

TYPED_TEST(journaled_sequence_group_fixture, constructibility)
{
    using type = typename TestFixture::type;

    ASSERT_TRUE(std::is_default_constructible_v<type>);
    ASSERT_TRUE(std::is_copy_constructible_v<type>);
    ASSERT_TRUE(std::is_nothrow_move_constructible_v<type>);
    ASSERT_TRUE(std::is_copy_assignable_v<type>);
    ASSERT_TRUE(std::is_nothrow_move_assignable_v<type>);
    ASSERT_TRUE(std::is_destructible_v<type>);
}

TYPED_TEST(journaled_sequence_group_fixture, construction_from_source)
{
    using type = typename TestFixture::type;

    EXPECT_NO_THROW((type{this->source(), 10}));
}

TYPED_TEST(journaled_sequence_group_fixture, construction_from_source_and_store)
{
    using type = typename TestFixture::type;

    EXPECT_NO_THROW((type{this->source(), this->generate_variants()}));
    EXPECT_NO_THROW((libjst::journaled_sequence_group{this->source(), this->generate_variants()}));
}

TYPED_TEST(journaled_sequence_group_fixture, clear)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, reset)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, at)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, size)
{
    using type = typename TestFixture::type;
    type js_group{};

    EXPECT_EQ(js_group.size(), 0u);
    js_group = type{this->source(), 10u};
    EXPECT_EQ(js_group.size(), 10u);
}

TYPED_TEST(journaled_sequence_group_fixture, load)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, save)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, begin)
{
    FAIL();
}

TYPED_TEST(journaled_sequence_group_fixture, end)
{
    FAIL();
}
