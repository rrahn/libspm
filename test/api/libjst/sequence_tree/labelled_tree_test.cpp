// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/sequence_tree/k_depth_tree.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>

#include "rcs_store_mock.hpp"

namespace jst::test::labelled_tree {

using source_t = std::string;
using variant_t = jst::test::variant<libjst::breakpoint, source_t, int, libjst::bit_vector<>>;

struct fixture {
    source_t source{};
    std::vector<variant_t> variants{};
    std::size_t coverage_size{};
    std::size_t window_size{};

    template <typename stream_t, typename this_t>
        requires std::same_as<std::remove_cvref_t<this_t>, fixture>
    friend stream_t & operator<<(stream_t & stream, this_t &&) {
        stream << "fixture";
        return stream;
    }
};

struct test : public ::testing::TestWithParam<fixture> {

    using rcs_store_t = mock_store<source_t>;

    rcs_store_t _mock;

    virtual void SetUp() override {
        _mock = rcs_store_t{GetParam().source, GetParam().coverage_size};
        std::ranges::for_each(GetParam().variants, [&] (auto var) {

            assert(std::ranges::size(libjst::coverage(var)) == _mock.size());
            _mock.insert(std::move(var));
        });
    }

    rcs_store_t const & get_mock() const noexcept{
        return _mock;
    }
};

} // namespace jst::test::labelled_tree

using namespace std::literals;


using node_id_type = std::pair<std::ptrdiff_t, bool>;

template <typename stream_t, typename id_t>
    requires std::same_as<std::remove_cvref_t<id_t>, node_id_type>
inline constexpr stream_t & operator<<(stream_t & stream, id_t && id) {
    stream << "<" << id.first << ", " << id.second << ">";
    return stream;
}

inline constexpr auto to_id = [] (std::ptrdiff_t left, std::ptrdiff_t right, bool is_alt) noexcept -> node_id_type {
    if (is_alt) {
        return node_id_type{left, true};
    } else {
        return node_id_type{left + right - 1, false};
    }
};

template <typename derived_t, typename extended_node_t>
class node_id_extension {
private:

    friend derived_t;

    node_id_extension() = default;


    constexpr void initialise() {
    };

    constexpr node_id_extension notify(extended_node_t const &) const {
        return {};
    }

    static constexpr derived_t const & as_derived(node_id_extension const & me) noexcept {
        return static_cast<derived_t const &>(me);
    }

public:

    constexpr node_id_type id() const noexcept {

        derived_t const & me = as_derived(*this);
        using seqan3::operator&;

        if ((as_derived(*this).get_second_breakpoint_id() & libjst::node_descriptor_id::second_first_right) != libjst::node_descriptor_id::nil) {
            return to_id(std::ranges::distance(me.rcs_store().variants().begin(), me.left_variant()),
                         std::ranges::distance(me.rcs_store().variants().begin(), me.left_variant()),
                         me.is_alt_node());
        } else {
            return to_id(std::ranges::distance(me.rcs_store().variants().begin(), me.left_variant()),
                         std::ranges::distance(me.rcs_store().variants().begin(), me.right_variant()),
                         me.is_alt_node());
        }
    }
};

template <typename base_tree_t>
using id_tree = libjst::extendable_tree<base_tree_t, node_id_extension>;

using fixture = jst::test::labelled_tree::fixture;
using variant_t = jst::test::labelled_tree::variant_t;
struct labelled_tree_test : public jst::test::labelled_tree::test
{
    using base_test_t = jst::test::labelled_tree::test;

    using typename base_test_t::rcs_store_t;
    using base_test_t::get_mock;
    using base_test_t::GetParam;

    using label_map_t = std::multimap<node_id_type, jst::test::labelled_tree::source_t>;

    label_map_t _label_map{};

    virtual void SetUp() override {
        base_test_t::SetUp();
        auto const & mock = get_mock();

        auto variant_beg = mock.variants().begin();
        auto variant_sent = mock.variants().end();
        auto const & src = mock.source();

        std::size_t prev_bp{};
        std::ptrdiff_t ref_id{0};
        std::ptrdiff_t var_id{1};
        for (auto variant_it = variant_beg; variant_it != variant_sent; ++variant_it, ++var_id, ++ref_id) {
            std::size_t nxt_bp = libjst::left_breakpoint(*variant_it);
            std::size_t prefix_span = nxt_bp - prev_bp;
            std::size_t var_span = libjst::breakpoint_span(*variant_it);

            _label_map.emplace_hint(_label_map.end(), node_id_type{ref_id, false}, src.substr(prev_bp, prefix_span));
            _label_map.emplace_hint(_label_map.end(), node_id_type{++ref_id, false}, src.substr(nxt_bp, var_span));
            _label_map.emplace_hint(_label_map.end(), node_id_type{var_id, true}, libjst::alt_sequence(*variant_it));
            prev_bp = nxt_bp + var_span;
        }
        _label_map.emplace_hint(_label_map.end(), to_id(var_id - 1, var_id, false), src.substr(prev_bp));
        std::cout << "Prebuild map:\n";
        print_map();
        std::cout << "\n";
    }

    void print_map() const noexcept {
        std::ranges::for_each(_label_map, [] (auto const & elem) {
            std::cout << elem.first << ": \"" << elem.second << "\"\n";
        });
    }

    template <typename tree_adator_t>
    auto make_tree(tree_adator_t && tree_adaptor) const noexcept {
        auto const & rcs_mock = get_mock();
        auto mock_tree = libjst::volatile_tree{rcs_mock} | libjst::labelled();
        using id_tree_t = id_tree<decltype(mock_tree)>;
        return tree_adaptor(id_tree_t{std::move(mock_tree)});
    }

    jst::test::labelled_tree::source_t const & expected_label_for(node_id_type const node_id) const noexcept {
        auto lbl_it = _label_map.find(node_id);
        EXPECT_TRUE(lbl_it != _label_map.end()) << "node_id = " << node_id;
        return lbl_it->second;
    }
};

// ----------------------------------------------------------------------------
// Test case definitions
// ----------------------------------------------------------------------------

TEST_P(labelled_tree_test, subtree_depth_0) {
    auto tree = make_tree(libjst::k_depth(0u));

    libjst::tree_traverser_base node_range{tree};
    for (auto && v : node_range) {
        std::cout << v.id() << ": " << "\"" << std::string{v.label().begin(), v.label().end()} << "\"\n";
        // EXPECT_RANGE_EQ(v.label(), expected_label_for(v.id()));
    }
}

TEST_P(labelled_tree_test, subtree_depth_1) {
    auto tree = make_tree(libjst::k_depth(1u));
    libjst::tree_traverser_base node_range{tree};
    for (auto && v : node_range) {
        std::cout << v.id() << ": " << "\"" << std::string{v.label().begin(), v.label().end()} << "\"\n";
        // EXPECT_RANGE_EQ(v.label(), expected_label_for(v.id()));
    }
}

TEST_P(labelled_tree_test, subtree_depth_2) {
    auto tree = make_tree(libjst::k_depth(2u));
    libjst::tree_traverser_base node_range{tree};
    for (auto && v : node_range) {
        std::cout << v.id() << ": " << "\"" << std::string{v.label().begin(), v.label().end()} << "\"\n";
        // EXPECT_RANGE_EQ(v.label(), expected_label_for(v.id()));
    }
}

// ----------------------------------------------------------------------------
// Test values
// ----------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(no_variant, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{},
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_first_base, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{0}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_last_base, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{7}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(snv_middle, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{4}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(two_snvs_scattered, labelled_tree_test, testing::Values(fixture{
           //  x  y
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{2}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}},
        variant_t{.position{5}, .insertion{"y"s}, .deletion{1}, .coverage{0, 1, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(two_snvs_next_to_each_other, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{3}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}},
        variant_t{.position{4}, .insertion{"y"s}, .deletion{1}, .coverage{0, 1, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));

INSTANTIATE_TEST_SUITE_P(two_snvs_at_same_breakpoint, labelled_tree_test, testing::Values(fixture{
    .source{"aaaabbbb"},
    .variants{
        variant_t{.position{4}, .insertion{"x"s}, .deletion{1}, .coverage{1, 0, 0, 0}},
        variant_t{.position{4}, .insertion{"y"s}, .deletion{1}, .coverage{0, 1, 0, 0}}
    },
    .coverage_size{4},
    .window_size{4}
}));
