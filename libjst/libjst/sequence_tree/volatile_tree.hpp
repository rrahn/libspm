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
#include <libjst/sequence_tree/empty_label.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>
#include <libjst/sequence_tree/rcs_store_node_base.hpp>
#include <libjst/referentially_compressed_sequence_store/rooted_rcs_store.hpp>

namespace libjst
{
    template <typename rcs_store_t>
    class volatile_tree {
    private:

        using rooted_rcs_store_type = rooted_rcs_store<rcs_store_t>;
        using variants_type = typename rooted_rcs_store_type::variant_map_type;
        using variant_type = std::ranges::range_value_t<variants_type>;

        class node_impl;

        rooted_rcs_store_type _rooted_rcs_store{};
        variant_type _bound{};

    public:

        volatile_tree() = default;
        volatile_tree(rcs_store_t const & rcs_store) noexcept :
            _rooted_rcs_store{rcs_store},
            _bound{*_rooted_rcs_store.variants().begin()}
        {
            using value_t = typename breakpoint::value_type;
            libjst::position(_bound) = breakpoint{static_cast<value_t>(std::ranges::size(_rooted_rcs_store.source()))};
        }

        constexpr node_impl root() const noexcept {
            auto first = std::ranges::begin(_rooted_rcs_store.variants());
            return node_impl{std::addressof(_rooted_rcs_store), first, std::ranges::next(first), _bound};
        }

        constexpr node_impl sink() const noexcept {
            auto sent = std::ranges::end(_rooted_rcs_store.variants());
            return node_impl{std::addressof(_rooted_rcs_store), sent, sent, _bound};
        }

        constexpr rooted_rcs_store_type const & data() const noexcept {
            return _rooted_rcs_store;
        }
    };

    template <typename rcs_store_t>
    class volatile_tree<rcs_store_t>::node_impl : public rcs_store_node_base<node_impl, rooted_rcs_store_type>  {
    private:

        friend volatile_tree;

        using base_t = rcs_store_node_base<node_impl, rooted_rcs_store_type>;
        using variant_iterator = typename rcs_node_traits<base_t>::variant_iterator;
        using variant_reference = std::iter_reference_t<variant_iterator>;

    private:

        variant_type const * _bound{};
        // constructor for the root/sink condition
        explicit constexpr node_impl(rooted_rcs_store_type const * rcs_store,
                                     variant_iterator left,
                                     variant_iterator right,
                                     variant_type const & bound) noexcept :
                base_t{rcs_store, std::move(left), std::move(right)},
                _bound{std::addressof(bound)}
        {
            assert(rcs_store != nullptr);
        }

    public:

        node_impl() = default;
        node_impl(node_impl const &) = default;
        node_impl(node_impl &&) = default;
        node_impl & operator=(node_impl const &) = default;
        node_impl & operator=(node_impl &&) = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return base_t::visit_next_alt();
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return base_t::visit_next_ref();
        }

        constexpr empty_label operator*() const noexcept {
            return {};
        }

    protected:
        constexpr variant_reference left_variant() const noexcept {
            return to_variant(base_t::get_left());
        }

        constexpr variant_reference right_variant() const noexcept {
            return to_variant(base_t::get_right());
        }

    private:

        constexpr variant_reference to_variant(variant_iterator const & it) const noexcept {
            if (it == base_t::sink())
                return *_bound;
            else
                return *it;
        }

        constexpr friend bool operator==(node_impl const & lhs, node_impl const & rhs) noexcept = default;
        // {
        //     return ;
        // }
    };
}  // namespace libjst
