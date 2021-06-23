// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <concepts>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>

#include <libjst/detail/journal_sequence_tree_event_branch.hpp>
#include <libjst/detail/journal_sequence_tree_event_join.hpp>
#include <libjst/detail/delta_event_shared.hpp>

using namespace std::literals;

namespace libjst::detail::test
{

using shared_delta_event_t = libjst::detail::delta_event_shared<char>;
using substitution_t = typename shared_delta_event_t::substitution_type;
using insertion_t = typename shared_delta_event_t::insertion_type;
using deletion_t = typename shared_delta_event_t::deletion_type;
using coverage_t = typename shared_delta_event_t::coverage_type;
using position_t = typename shared_delta_event_t::position_type;

enum struct event_category
{
    branch,
    join
};

} // namespace libjst::detail::test

struct jst_event_fixture
{
    libjst::detail::test::shared_delta_event_t event{};
    libjst::detail::test::position_t expected_position{};
    libjst::detail::test::event_category category{};
};

struct jst_event_test : public ::testing::TestWithParam<jst_event_fixture>
{
    using shared_delta_event_t = libjst::detail::test::shared_delta_event_t;
    using branch_event_t = libjst::detail::journal_sequence_tree_event_branch<shared_delta_event_t>;
    using join_event_t = libjst::detail::journal_sequence_tree_event_join<shared_delta_event_t>;
    using jst_event_variant_t = std::variant<branch_event_t, join_event_t>;

    jst_event_variant_t test_event{};
    shared_delta_event_t expected_event{};
    libjst::detail::test::coverage_t expected_coverage;
    libjst::detail::test::position_t expected_position{};

    void SetUp() override
    {
        expected_event = GetParam().event;
        if (GetParam().category == libjst::detail::test::event_category::branch)
            test_event = branch_event_t{std::addressof(expected_event)};
        else
            test_event = join_event_t{std::addressof(expected_event)};
        expected_position = GetParam().expected_position;
        expected_coverage = GetParam().event.coverage();
    }

    template <typename event_variant_t, typename fn_t>
    decltype(auto) apply(event_variant_t && event_variant, fn_t && fn) const
    {
        return std::visit([&] (auto && event) { return fn(event); }, std::forward<event_variant_t>(event_variant));
    }

    jst_event_variant_t get_testable_event(shared_delta_event_t * event_ptr) const noexcept
    {
        if (std::holds_alternative<branch_event_t>(test_event))
            return branch_event_t{event_ptr};
        else
            return join_event_t{event_ptr};
    }
};

TEST_P(jst_event_test, construction)
{
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<branch_event_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<branch_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<branch_event_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<branch_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<branch_event_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<branch_event_t>);
    EXPECT_TRUE((std::is_nothrow_constructible_v<branch_event_t, shared_delta_event_t *>));

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<join_event_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<join_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<join_event_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<join_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<join_event_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<join_event_t>);
    EXPECT_TRUE((std::is_nothrow_constructible_v<join_event_t, shared_delta_event_t *>));
}

TEST_P(jst_event_test, position)
{
    EXPECT_EQ((apply(test_event, [] (auto & e) { return e.position();})), expected_position);
    EXPECT_EQ((apply(std::as_const(test_event), [] (auto const & e) { return e.position();})), expected_position);
}

TEST_P(jst_event_test, coverage)
{
    EXPECT_EQ((apply(test_event, [] (auto & e) { return e.coverage();})), expected_coverage);
    EXPECT_EQ((apply(std::as_const(test_event), [] (auto const & e) { return e.coverage();})), expected_coverage);
}

TEST_P(jst_event_test, event_handle)
{
    EXPECT_EQ(*(apply(test_event, [] (auto & e) { return e.event_handle();})), expected_event);
    EXPECT_EQ(*(apply(std::as_const(test_event), [] (auto const & e) { return e.event_handle();})), expected_event);
}

TEST_P(jst_event_test, equality)
{
    libjst::detail::test::shared_delta_event_t other_delta_event{expected_event};
    jst_event_variant_t event_pointing_to_same_address = get_testable_event(std::addressof(expected_event));
    jst_event_variant_t event_pointing_to_another_address = get_testable_event(std::addressof(other_delta_event));

    EXPECT_TRUE(std::equality_comparable<branch_event_t>);
    EXPECT_TRUE(std::equality_comparable<join_event_t>);
    EXPECT_EQ(test_event, test_event);
    EXPECT_EQ(test_event, event_pointing_to_same_address);
    EXPECT_NE(test_event, event_pointing_to_another_address);
}

TEST_P(jst_event_test, ordering_by_different_event_position)
{
    // Create an event that is greater than the given one but assume that the event position is less than the maximal
    // number of size_t.
    ASSERT_LT(expected_position.offset, std::numeric_limits<size_t>::max());
    libjst::detail::test::shared_delta_event_t greater_delta_event{expected_position + 1,
                                                                   expected_event.delta_variant(),
                                                                   expected_event.coverage()};
    jst_event_variant_t event_pointing_to_same_event = get_testable_event(std::addressof(expected_event));
    jst_event_variant_t event_pointing_to_greater_event = get_testable_event(std::addressof(greater_delta_event));

    { // less and less equal
        EXPECT_LT(std::as_const(test_event), std::as_const(event_pointing_to_greater_event));
        EXPECT_LE(std::as_const(test_event), std::as_const(event_pointing_to_same_event));
        EXPECT_LE(std::as_const(test_event), std::as_const(test_event));
    }

    { // greater and greater equal
        EXPECT_GT(event_pointing_to_greater_event, test_event);
        EXPECT_GE(event_pointing_to_same_event, test_event);
        EXPECT_GE(test_event, test_event);
    }
}

TEST_P(jst_event_test, ordering_by_same_event_position_with_different_delta_kind)
{
    using namespace std::literals;

    // Create events of different delta kinds at the same position and compare ordering of jst events based
    // on the delta kind, i.e. insertion < substitution < deletion for branch events and
    // deletion < substitution < insertion for join events.
    // The following events for comparison are adapted such that they have the same position as the test event.
    // Note that the join event position is position of the delta event + deletion size.
    bool const is_join_event = test_event.index() == 1;
    libjst::detail::test::shared_delta_event_t insertion_event{expected_position,
                                                               libjst::detail::test::insertion_t{"a"s},
                                                               expected_event.coverage()};
    jst_event_variant_t event_to_insertion = get_testable_event(std::addressof(insertion_event));

    libjst::detail::test::shared_delta_event_t deletion_event{expected_position - (is_join_event ? 1 : 0),
                                                              libjst::detail::test::deletion_t{1},
                                                              expected_event.coverage()};
    jst_event_variant_t event_to_deletion = get_testable_event(std::addressof(deletion_event));

    libjst::detail::test::shared_delta_event_t substitution_event{expected_position - (is_join_event ? 1 : 0),
                                                                  libjst::detail::test::substitution_t{"a"s},
                                                                  expected_event.coverage()};
    jst_event_variant_t event_to_substitution = get_testable_event(std::addressof(substitution_event));

    // The set index of the delta variant of the test event: 0 = insertion, 2 = substitution, 3 = deletion
    int32_t const test_index = apply(test_event, [] (auto & e) { return e.event_handle(); })->delta_variant().index();

    // Returns the result of invoking either cmp_branch or cmp_join with the test index and the expected one.
    auto with_index = [&] (auto cmp_branch, auto cmp_join, int32_t expected_index)
    {
        return (is_join_event) ? cmp_join(test_index, expected_index) : cmp_branch(test_index, expected_index);
    };

    { // less and less equal
        EXPECT_EQ(std::as_const(test_event) < event_to_insertion,
                  with_index(std::less<int32_t>{}, std::greater<int32_t>{}, 0));
        EXPECT_EQ(std::as_const(test_event) < event_to_substitution,
                  with_index(std::less<int32_t>{}, std::greater<int32_t>{}, 2));
        EXPECT_EQ(std::as_const(test_event) < event_to_deletion,
                  with_index(std::less<int32_t>{}, std::greater<int32_t>{}, 3));

        EXPECT_EQ(std::as_const(test_event) <= event_to_insertion,
                  with_index(std::less_equal<int32_t>{}, std::greater_equal<int32_t>{}, 0));
        EXPECT_EQ(std::as_const(test_event) <= event_to_substitution,
                  with_index(std::less_equal<int32_t>{}, std::greater_equal<int32_t>{}, 2));
        EXPECT_EQ(std::as_const(test_event) <= event_to_deletion,
                  with_index(std::less_equal<int32_t>{}, std::greater_equal<int32_t>{}, 3));
    }

    { // greater and greater equal
        EXPECT_EQ(std::as_const(test_event) > event_to_insertion,
                  with_index(std::greater<int32_t>{}, std::less<int32_t>{}, 0));
        EXPECT_EQ(std::as_const(test_event) > event_to_substitution,
                  with_index(std::greater<int32_t>{}, std::less<int32_t>{}, 2));
        EXPECT_EQ(std::as_const(test_event) > event_to_deletion,
                  with_index(std::greater<int32_t>{}, std::less<int32_t>{}, 3));

        EXPECT_EQ(std::as_const(test_event) >= event_to_insertion,
                  with_index(std::greater_equal<int32_t>{}, std::less_equal<int32_t>{}, 0));
        EXPECT_EQ(std::as_const(test_event) >= event_to_substitution,
                  with_index(std::greater_equal<int32_t>{}, std::less_equal<int32_t>{}, 2));
        EXPECT_EQ(std::as_const(test_event) >= event_to_deletion,
                  with_index(std::greater_equal<int32_t>{}, std::less_equal<int32_t>{}, 3));
    }
}

TEST_P(jst_event_test, ordering_by_position)
{
    // Get position that is greater than the one from the test event.
    ASSERT_LT(expected_position.offset, std::numeric_limits<size_t>::max());
    libjst::detail::test::position_t const greater_position = expected_position + 1;

    { // less and less than equal
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return e < greater_position; }));
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return e <= greater_position; }));
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return e <= expected_position; }));
    }

    { // greater and greater than equal
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return greater_position > e; }));
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return greater_position >= e; }));
        EXPECT_TRUE(apply(std::as_const(test_event), [&] (auto & e) { return greater_position >= e; }));
    }
}
