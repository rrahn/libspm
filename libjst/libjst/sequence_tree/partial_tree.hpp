// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides merge tree with maximal branch-free nodes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/breakend_site_trimmed.hpp>
#include <libjst/sequence_tree/breakend_site.hpp>
#include <libjst/sequence_tree/breakpoint_node.hpp>
#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/empty_label.hpp>
#include <libjst/variant/breakpoint.hpp>

namespace libjst
{
    template <typename rcs_store_t>
    class partial_tree {
    private:

        using variants_type = typename rcs_store_t::variant_map_type;
        using variant_type = std::ranges::range_value_t<variants_type>;
        using breakend_iterator = std::ranges::iterator_t<variants_type const &>;
        using base_node_type = breakpoint_node<breakend_iterator>;
        using position_type = typename base_node_type::position_type;
        using trimmed_position_type = breakend_site_trimmed<position_type>;
        using position_value_type = typename trimmed_position_type::position_value_type;

        class node_impl;

        rcs_store_t const & _rcs_store{};
        position_type _low_base{};
        trimmed_position_type _partial_low_nil{};
        trimmed_position_type _partial_high_nil{};

    public:

        using size_type = position_value_type;
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr partial_tree() = default; //!< Default.

        partial_tree(rcs_store_t const & rcs_store, position_value_type root_position, position_value_type count) noexcept :
            _rcs_store{rcs_store}
        {
            // initiate low boundary
            auto low = std::ranges::lower_bound(std::ranges::next(std::ranges::begin(data().variants())),
                                                std::ranges::end(data().variants()),
                                                root_position,
                                                std::ranges::less{},
                                                [&] (auto breakend_proxy) -> position_value_type {
                                                     return libjst::position(std::move(breakend_proxy));
                                                });
            assert(root_position <= static_cast<position_value_type>(libjst::position(*low)));
            _low_base = position_type{std::ranges::prev(low), breakpoint_end::low};

            auto high_base_it = std::ranges::prev(std::ranges::end(data().variants()));
            _partial_low_nil = trimmed_position_type{position_type{high_base_it, breakpoint_end::low}, root_position};
            _partial_high_nil = trimmed_position_type{position_type{high_base_it, breakpoint_end::high},
                                                      root_position + count};
        }
        //!\}

        constexpr node_impl root() const noexcept {
            base_node_type base_root{_low_base, _low_base};
            return node_impl{base_root.next_ref(), _partial_low_nil, _partial_high_nil};
        }

        constexpr nil_node_t sink() const noexcept {
            return nil_node;
        }

        constexpr rcs_store_t const & data() const noexcept {
            return _rcs_store;
        }
   };

    template <typename rcs_store_t, std::integral offset_t, std::integral count_t>
    partial_tree(rcs_store_t const &, offset_t, count_t) -> partial_tree<rcs_store_t>;

    template <typename base_tree_t>
    class partial_tree<base_tree_t>::node_impl : public base_node_type {
    private:

        friend partial_tree;

        trimmed_position_type _partial_lowest;
        trimmed_position_type _partial_highest;

        explicit constexpr node_impl(base_node_type base_node,
                                     trimmed_position_type partial_lowest,
                                     trimmed_position_type partial_highest) noexcept :
            base_node_type{std::move(base_node)},
            _partial_lowest{std::move(partial_lowest)},
            _partial_highest{std::move(partial_highest)}
        {
        }

    public:

        using low_position_type = trimmed_position_type;
        using high_position_type = trimmed_position_type;

        constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (reached_highest() /*&& (libjst::position(_partial_highest) != libjst::position(_partial_highest.base()))*/) {
                base_node_type child{base_node_type::low_boundary(), base_node_type::high_boundary()};
                child.toggle_alternate_path();
                return node_impl{std::move(child), _partial_highest, trimmed_position_type{_partial_highest.base()}};
            }

            if (auto maybe_child = base_node_type::next_alt(); maybe_child.has_value()) {
                trimmed_position_type child_highest = (!this->on_alternate_path()) ?
                                                      trimmed_position_type{_partial_highest.base()} :
                                                      _partial_highest;
                return node_impl{std::move(*maybe_child), _partial_lowest, std::move(child_highest)};
            }
            return std::nullopt;
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (is_leaf()) {
                return std::nullopt;
            }
            return node_impl{base_node_type::next_ref(), _partial_lowest, _partial_highest};
        }

        constexpr empty_label operator*() const noexcept {
            return {};
        }

        constexpr low_position_type low_boundary() const noexcept {
            position_type low_base = base_node_type::low_boundary();
            if (libjst::position(low_base) < libjst::position(_partial_lowest)) {
                return _partial_lowest;
            }
            return low_position_type{base_node_type::low_boundary()};
        }

        constexpr high_position_type high_boundary() const noexcept {
            position_type high_base = base_node_type::high_boundary();
            if (reached_highest()) {
                return _partial_highest;
            }
            return low_position_type{std::move(high_base)};
        }

        constexpr bool is_leaf() const noexcept {
            return reached_highest() || (base_node_type::high_boundary() == _partial_highest.base());
        }
    private:

        constexpr bool reached_highest() const noexcept {
            return !this->on_alternate_path() &&
                   libjst::position(_partial_highest) <= libjst::position(base_node_type::high_boundary());
        }

        constexpr friend bool operator==(node_impl const & lhs, nil_node_t const &) noexcept {
            return lhs.is_leaf();
        }
    };
}  // namespace libjst
