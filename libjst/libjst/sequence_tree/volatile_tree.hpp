// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/empty_label.hpp>
#include <libjst/sequence_tree/breakpoint_node.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>
#include <libjst/sequence_tree/rcs_store_node_base.hpp>

namespace libjst
{
    template <typename rcs_store_t>
    class volatile_tree {
    private:

        using rooted_rcs_store_type = rcs_store_t;
        using variants_type = typename rooted_rcs_store_type::variant_map_type;
        using variant_type = std::ranges::range_value_t<variants_type>;
        using breakend_iterator = std::ranges::iterator_t<variants_type const &>;
        using node_type = breakpoint_node<breakend_iterator>;
        using position_type = typename node_type::position_type;

        class node_impl;

        rooted_rcs_store_type const & _rooted_rcs_store{};
        position_type _low_nil{};
        position_type _high_nil{};

    public:

        volatile_tree() = delete;
        volatile_tree(rcs_store_t const & rcs_store) noexcept :
            _rooted_rcs_store{rcs_store}
        {
            _low_nil = position_type{std::ranges::begin(data().variants()), breakpoint_end::low};
            _high_nil = position_type{std::ranges::prev(std::ranges::end(data().variants())), breakpoint_end::high};
        }

        constexpr node_impl root() const noexcept {
            node_type root_base{_low_nil, _low_nil};
            return node_impl{root_base.next_ref(), _high_nil};
        }

        constexpr nil_node_t sink() const noexcept {
            return nil_node;
        }

        constexpr rooted_rcs_store_type const & data() const noexcept {
            return _rooted_rcs_store;
        }
    };

    template <typename rcs_store_t>
    class volatile_tree<rcs_store_t>::node_impl : public node_type {
    private:

        friend volatile_tree;

        using base_t = node_type;

    private:

        position_type _nil{};

        explicit constexpr node_impl(base_t && base_node, position_type nil) noexcept :
                base_t{std::move(base_node)},
                _nil{std::move(nil)}
        {
        }

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (auto child = base_t::next_alt(); child.has_value())
                return node_impl{std::move(*child), _nil};
            return std::nullopt;
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (is_leaf()) // never reaching this?
                return std::nullopt;
            return node_impl{base_t::next_ref(), _nil};
        }

        constexpr empty_label operator*() const noexcept {
            return {};
        }

    private:

        constexpr bool is_leaf() const noexcept {
            return base_t::high_boundary() == _nil;
        }

        constexpr friend bool operator==(node_impl const & lhs, nil_node_t const &) noexcept
        {
            return lhs.is_leaf();
        }
    };

    namespace _tree_factory {
        inline constexpr struct _make_volatile
        {
            template <typename rcs_store_t>
            constexpr auto operator()(rcs_store_t const & rcs_store) const
                noexcept(std::is_nothrow_constructible_v<volatile_tree<rcs_store_t>>)
                -> volatile_tree<rcs_store_t>
            {
                return volatile_tree<rcs_store_t>{rcs_store};
            }

            template <typename ...args_t>
                requires (sizeof...(args_t) == 0)
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _make_volatile, args_t...>)
                -> jst::contrib::closure_result_t<_make_volatile, args_t...>
            {
                return jst::contrib::make_closure(_make_volatile{}, (args_t&&) args...);
            }
        } make_volatile{};
    } // namespace _tree_factory

    using _tree_factory::make_volatile;

}  // namespace libjst
