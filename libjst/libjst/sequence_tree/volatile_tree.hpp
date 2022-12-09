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

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>
#include <libjst/sequence_tree/rcs_store_node_base.hpp>
#include <libjst/referentially_compressed_sequence_store/rooted_rcs_store.hpp>

namespace libjst
{
    template <typename rcs_store_t>
    class volatile_tree {
    private:

        using rooted_rcs_store_type = rooted_rcs_store<rcs_store_t>;

        class node_impl;

        rooted_rcs_store_type _rooted_rcs_store{};
    public:

        volatile_tree(rcs_store_t const & rcs_store) noexcept : _rooted_rcs_store{rcs_store}
        {}

        constexpr node_impl root() const noexcept {
            auto first = std::ranges::begin(_rooted_rcs_store.variants());
            return node_impl{_rooted_rcs_store, first, std::ranges::next(first)};
        }

        constexpr node_impl sink() const noexcept {
            auto sent = std::ranges::end(_rooted_rcs_store.variants());
            return node_impl{_rooted_rcs_store, sent, sent};
        }
    };

    template <typename rcs_store_t>
    class volatile_tree<rcs_store_t>::node_impl : public rcs_store_node_base<node_impl, rooted_rcs_store_type>  {
    private:

        friend volatile_tree;

        using base_t = rcs_store_node_base<node_impl, rooted_rcs_store_type>;
        using variant_iterator = typename rcs_node_traits<base_t>::variant_iterator;

    private:
        // constructor for the root/sink condition
        explicit constexpr node_impl(rooted_rcs_store_type const & rcs_store,
                                     variant_iterator left,
                                     variant_iterator right) noexcept :
                base_t{rcs_store, std::move(left), std::move(right)}
        {}

    public:

        explicit constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return base_t::visit_next_alt();
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return base_t::visit_next_ref();
        }

    private:

        constexpr friend bool operator==(node_impl const & lhs, node_impl const & rhs) noexcept {
            return lhs.left_variant() == rhs.left_variant() && lhs.right_variant() == rhs.right_variant();
        }
    };
}  // namespace libjst
