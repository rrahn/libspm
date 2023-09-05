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

#include <libjst/sequence_tree/breakend_site_min.hpp>
#include <libjst/sequence_tree/breakend_site_partial.hpp>
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
        using breakend_type = decltype(std::declval<position_type const &>().get_breakend());
        using partial_position_type = breakend_site_partial<breakend_type>;
        using _low_position_type = breakend_site_min<partial_position_type>;
        using _high_position_type = breakend_site_trimmed<partial_position_type>;
        using position_value_type = typename _high_position_type::position_value_type;

        class node_impl;

        rcs_store_t const & _rcs_store;
        position_type _low_base{};
        _low_position_type _partial_low_nil{};
        _high_position_type _partial_high_nil{};
    protected:

        constexpr explicit partial_tree(rcs_store_t const & rcs_store) noexcept :
            _rcs_store{rcs_store}
        {
        }

    public:

        using size_type = position_value_type;
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr partial_tree() = default; //!< Default.

        partial_tree(rcs_store_t const & rcs_store, position_value_type root_position, position_value_type count) noexcept :
            partial_tree{rcs_store}
        {
            // initiate low boundary
            position_value_type end_position = std::min<position_value_type>(root_position + count,
                                                                             rcs_store.source().size());
            auto low = std::ranges::lower_bound(std::ranges::next(std::ranges::begin(data().variants())),
                                                std::ranges::end(data().variants()),
                                                root_position,
                                                std::ranges::less{},
                                                [&] (auto breakend_proxy) -> position_value_type {
                                                     return libjst::position(std::move(breakend_proxy));
                                                });

            assert(root_position <= static_cast<position_value_type>(libjst::position(*low)));

            set_low_base(position_type{std::ranges::prev(low), breakpoint_end::low});
            auto high_base_it = std::ranges::prev(std::ranges::end(data().variants()));
            // set low nil node:
            partial_position_type partial_root{_low_base.get_breakend(), high_base_it, breakpoint_end::low};
            set_low_nil(_low_position_type{std::move(partial_root), root_position});

            auto high = std::ranges::lower_bound(low,
                                                 high_base_it,
                                                 end_position,
                                                 std::ranges::less{},
                                                 [&] (auto breakend_proxy) -> position_value_type {
                                                      return libjst::position(std::move(breakend_proxy));
                                                 });

            assert(end_position <= static_cast<position_value_type>(libjst::position(*high)));
            partial_position_type partial_sink{high, high_base_it, breakpoint_end::high};
            set_high_nil(_high_position_type{std::move(partial_sink), end_position});
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
    protected:

        constexpr void set_low_base(position_type low_base) noexcept {
            _low_base = std::move(low_base);
        }

        constexpr void set_low_nil(_low_position_type low_nil) noexcept {
            _partial_low_nil = std::move(low_nil);
        }

        constexpr void set_high_nil(_high_position_type high_nil) noexcept {
            _partial_high_nil = std::move(high_nil);
        }
    };

    template <typename rcs_store_t, std::integral offset_t, std::integral count_t>
    partial_tree(rcs_store_t const &, offset_t, count_t) -> partial_tree<rcs_store_t>;

    template <typename base_tree_t>
    class partial_tree<base_tree_t>::node_impl : public base_node_type {
    private:

        friend partial_tree;

        _low_position_type _partial_lowest;
        _high_position_type _partial_highest;
        bool _passed_high_bound{false};

        explicit constexpr node_impl(base_node_type base_node,
                                     _low_position_type partial_lowest,
                                     _high_position_type partial_highest,
                                     bool passed_high_bound = false) noexcept :
            base_node_type{std::move(base_node)},
            _partial_lowest{std::move(partial_lowest)},
            _partial_highest{std::move(partial_highest)},
            _passed_high_bound{passed_high_bound}
        {
        }

    public:

        using low_position_type = _low_position_type;
        using high_position_type = _high_position_type;

        constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (reached_highest() && !_passed_high_bound) {
                assert(libjst::position(base_node_type::low_boundary()) <= libjst::position(_partial_highest));
                assert(libjst::position(_partial_highest) <= libjst::position(base_node_type::high_boundary()));

                base_node_type child{base_node_type::low_boundary(), base_node_type::high_boundary()};
                child.toggle_alternate_path();
                partial_position_type _low{base_node_type::low_boundary().get_breakend(),
                                           _partial_lowest.base().get_bound(),
                                           base_node_type::low_boundary().get_breakend_site()};
                partial_position_type global_bound{_partial_highest.base().get_bound(),
                                                   _partial_highest.base().get_bound(),
                                                   _partial_highest.get_breakend_site()};
                return node_impl{std::move(child),
                                 low_position_type{std::move(_low), libjst::position(_partial_highest)},
                                 high_position_type{std::move(global_bound)},
                                 true};
            }

            if (auto maybe_child = base_node_type::next_alt(); maybe_child.has_value()) {
                high_position_type child_highest{_partial_highest};
                if (!this->on_alternate_path()) {
                        child_highest = high_position_type{partial_position_type{_partial_highest.base().get_bound(),
                                                                                 _partial_highest.base().get_bound(),
                                                                                 _partial_highest.get_breakend_site()}};
                }
                return node_impl{std::move(*maybe_child), _partial_lowest, std::move(child_highest), _passed_high_bound};
            }
            return std::nullopt;
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (is_leaf()) {
                return std::nullopt;
            }
            return node_impl{base_node_type::next_ref(), _partial_lowest, _partial_highest, _passed_high_bound};
        }

        constexpr empty_label operator*() const noexcept {
            return {};
        }

        constexpr low_position_type low_boundary() const noexcept {
            position_type low_base = base_node_type::low_boundary();
            if (libjst::position(low_base) < libjst::position(_partial_lowest)) {
                return _partial_lowest;
            }
            return low_position_type{partial_position_type{std::move(low_base)}};
        }

        constexpr high_position_type high_boundary() const noexcept {
            position_type high_base = base_node_type::high_boundary();
            if (reached_highest()) {
                return _partial_highest;
            }
            return high_position_type{partial_position_type{std::move(high_base)}};
        }

    protected:

        constexpr bool is_leaf() const noexcept {
            auto high_bound = base_node_type::high_boundary();
            partial_position_type nil{_partial_highest.base().get_breakend(),
                                      high_bound.get_breakend(),
                                      high_bound.get_breakend_site()};
            return reached_highest() || (nil == _partial_highest.base());
        }


        template <typename breakend_site_t>
        constexpr void reset_low(breakend_site_t new_low) {
            base_node_type tmp{new_low, new_low};
            static_cast<base_node_type &>(*this) = tmp.next_ref();
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
