// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides base implementation for the sequence tree nodes based on the rcs_store.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <optional>
#include <ranges>

#include <libjst/variant/concept.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>
// #include <libjst/sequence_tree/breakpoint_node.hpp>
namespace libjst
{
    template <typename derived_t, typename rcs_store_t>
    class rcs_store_node_base : public node_descriptor {
    public:
        using rcs_store_type = rcs_store_t;
    protected:

        friend derived_t;

        using variant_map_type = typename rcs_store_t::variant_map_type;
        using variant_iterator = std::ranges::iterator_t<variant_map_type const>;
        using difference_type =  std::iter_difference_t<variant_iterator>;
        using variant_type = std::iter_value_t<variant_iterator>;
        using breakpoint_type = libjst::variant_breakpoint_t<variant_type>;
        using optional_node_type = std::optional<derived_t>;
        // using base_node_type = breakpoint_node<variant_iterator>;
        using position_type = typename libjst::variant_position_t<variant_type>::value_type;

    private:
        rcs_store_type const * _rcs_store{};
        variant_iterator _left_variant{};
        variant_iterator _right_variant{};
        variant_iterator _next_variant{};

        // base_node_type _base{};

    protected:

        rcs_store_node_base() = default;
        constexpr rcs_store_node_base(rcs_store_type const * rcs_store,
                                      variant_iterator left_variant,
                                      variant_iterator right_variant) noexcept :
            _rcs_store{rcs_store},
            _left_variant{std::move(left_variant)},
            _right_variant{std::move(right_variant)}
        {
            assert(_rcs_store != nullptr);
            set_next(next_variant_after(_right_variant));
            initialise_reference_state_from();
        }

        // only from branching node!
        constexpr optional_node_type visit_next_alt() const noexcept {
            if (is_ref_node() && is_branching()) {
                derived_t child{as_derived(*this)};
                child.set_left(get_right());
                child.activate_state(node_state::variant);
                return child;
            } else {
                return std::nullopt;
            }
        }

        constexpr optional_node_type visit_next_ref() const noexcept {
            derived_t child{as_derived(*this)};
            if (node_descriptor::from_reference()) {
                if (is_leaf_of_alternate_subtree()) {
                    return std::nullopt;
                } else {
                    return visit_next_ref_impl(std::move(child));
                }
            } else {
                assert(get_left() == get_right());
                child.set_left(get_right());
                child.set_right(find_next_valid_right_variant());
                child.set_next(next_variant_after(child.get_right()));
                child.initialise_reference_state_from();
            }
            return child;
        }

        constexpr bool is_ref_node() const noexcept {
            return node_descriptor::from_reference();
        }

        constexpr bool is_alt_node() const noexcept {
            return node_descriptor::from_variant();
        }

        constexpr bool on_alternate_path() const noexcept {
            return node_descriptor::on_alternate_path();
        }

        /*!\brief Returns the breakpoint of the left end of the current node.
         *
         * The correct position of the current node depends on the node type. If `this` is a branching node, i.e.
         * libjst::rcs_store_node_base::is_branching_node() evaluates to `true`, this member returns the end position
         * of the left variant, otherwise it returns the begin position of the left variant.
         *
         * \returns Source position of the left end of the current node.
         */
        constexpr position_type low_breakend() const {
            if (get_left() == sink())
                return std::ranges::size(rcs_store().source());

            auto bounded_right_position = [&] (variant_iterator const it) -> position_type {
                if (it == sink()) {
                    return std::ranges::size(rcs_store().source());
                } else {
                    return libjst::position(*it).value();
                }
            };

            return (node_descriptor::left_break().from_left_begin())
                        ? libjst::left_breakpoint(*get_left()).value()
                        : std::min<position_type>(libjst::right_breakpoint(*get_left()).value(), bounded_right_position(get_right()));
        }

        /*!\brief Returns the breakpoint of the right end of the current node.
         *
         * The correct position of the current node depends on the node type. If `this` is not a branching node, i.e.
         * libjst::rcs_store_node_base::is_branching_node() evaluates to `false`, this member returns the position of
         * the right variant, otherwise it returns the end position of the left variant.
         *
         * \returns Source position of the right end of the current node.
         */
        constexpr position_type high_breakend() const {

            auto bounded_right_breakpoint = [&] (variant_iterator const it) -> position_type {
                if (it == sink()) {
                    return std::ranges::size(rcs_store().source());
                } else {
                    return libjst::right_breakpoint(*it).value();
                }
            };

            auto bounded_left_breakpoint = [&] (variant_iterator const it) -> position_type {
                if (it == sink()) {
                    return std::ranges::size(rcs_store().source());
                } else {
                    return libjst::left_breakpoint(*it).value();
                }
            };

            using seqan3::operator&;
            if (node_descriptor::right_break().from_left_end()) {
                return bounded_right_breakpoint(get_left());
            } else if (node_descriptor::right_break().from_right_end()) {
                return bounded_right_breakpoint(get_right());
            } else {
                return bounded_left_breakpoint(get_right());
            }
        }

        constexpr rcs_store_t const & rcs_store() const noexcept {
            assert(_rcs_store != nullptr);
            return *_rcs_store;
        }

        constexpr void set_left(variant_iterator new_left) noexcept {
            _left_variant = std::move(new_left);
        }

        constexpr variant_iterator get_left() const noexcept {
            return _left_variant;
        }

        constexpr void set_right(variant_iterator new_right) noexcept {
            _right_variant = std::move(new_right);
        }

        constexpr variant_iterator get_right() const noexcept {
            return _right_variant;
        }

        // for now we only consider the SNVs!
        constexpr variant_iterator next_variant_after(variant_iterator it) const noexcept {
            return std::ranges::find_if(it, sink(), [it] (auto && variant) {
                return libjst::left_breakpoint(*it) < libjst::left_breakpoint(variant);
            });
        }

        constexpr void set_next(variant_iterator new_next) noexcept {
            _next_variant = std::move(new_next);
        }

        constexpr bool is_branching() const noexcept {
            return node_descriptor::is_branching();
        }

        constexpr variant_iterator sink() const noexcept {
            return std::ranges::prev(std::ranges::end(rcs_store().variants()));
        }

        constexpr bool is_nil() const noexcept {
            return from_reference() && get_right() == sink() && !node_descriptor::right_break().from_left_end();
        }

        constexpr void initialise_reference_state_from() noexcept {
            bool const is_branching = (get_right() != sink()) && libjst::position(*get_right()).is_left_end();
            bool const is_last = is_branching && right_before_next();

            if (is_branching) {
                if (is_last) {
                    node_descriptor::activate_state(node_state::last_branching_after_left_end);
                } else {
                    node_descriptor::activate_state(node_state::branching_after_left_end);
                }
            } else {
                node_descriptor::activate_state(node_state::non_branching_after_left);
            }
        }

    private:

        constexpr bool is_leaf_of_alternate_subtree() const noexcept {
            return node_descriptor::on_alternate_path() &&
                   get_right() == sink() &&
                   !node_descriptor::right_break().from_left_end();
        }

        constexpr variant_iterator get_next() const noexcept {
            return _next_variant;
        }

        constexpr variant_iterator find_next_valid_right_variant() const noexcept {
            variant_iterator next_right = get_right();
            assert(is_alt_node());
            auto const min_ref_position = libjst::right_breakpoint(*next_right);
            // TODO: Do we need the ordering here?
            next_right = std::ranges::find_if_not(std::ranges::next(next_right), sink(), [&] (auto && var) {
                return libjst::left_breakpoint(var) < min_ref_position;
            });
            assert(next_right != get_right());
            return next_right;
        }

        static constexpr derived_t visit_next_ref_impl(derived_t child) noexcept {
            node_state parent_node_state = static_cast<node_state>(child);
            // Move through the states.
            switch(parent_node_state) {
                case node_state::branching_after_left_end: [[fallthrough]];
                case node_state::branching_after_left_begin:
                    child.set_left(child.get_right());
                    child.set_right(std::ranges::next(child.get_right()));
                    break;
                case node_state::last_branching_after_left_end: [[fallthrough]];
                case node_state::last_branching_after_left_begin: [[fallthrough]];
                case node_state::non_branching_after_left: [[fallthrough]];
                case node_state::non_branching_including_left:
                    child.set_left(child.get_right());
                    child.set_right(child.get_next());
                    child.set_next(child.next_variant_after(child.get_next()));
                    break;
                default: /*no-op*/;
            }

            // condition: in last node parent can only be in G or H
            auto bounded_left_breakpoint = [&] (variant_iterator it) {
                if (it == child.sink()) {
                    return breakpoint{static_cast<uint32_t>(std::ranges::size(child.rcs_store().source()))};
                } else {
                    return libjst::left_breakpoint(*it);
                }
            };

            auto bounded_breakpoint = [&] (variant_iterator const it) {
                if (it == child.sink()) {
                    return breakpoint{static_cast<uint32_t>(std::ranges::size(child.rcs_store().source()))};
                } else {
                    return libjst::position(*it);
                }
            };

            // Update the breakpoint ids of the child node.
            switch(parent_node_state) {
                  case node_state::last_branching_after_left_end: { [[fallthrough]];
                } case node_state::last_branching_after_left_begin: {
                    if (libjst::right_breakpoint(*child.get_left()) < bounded_left_breakpoint(child.get_right())) { // => {C/D}
                        if (child.get_right() != child.sink() && libjst::position(*child.get_right()).is_left_end()) {
                            child.activate_state(node_state::last_non_branching_left_only);
                        } else {
                            child.activate_state(node_state::non_branching_left_only);
                        }
                    } else if (bounded_breakpoint(child.get_right()).is_left_end()) { // => {E/F}
                        if (child.right_before_next()) {
                            child.activate_state(node_state::last_branching_after_left_begin);
                        } else {
                            child.activate_state(node_state::branching_after_left_begin);
                        }
                    } else { // => {H}
                        child.activate_state(node_state::non_branching_including_left);
                    }
                    break;
                } case node_state::last_non_branching_left_only: { // C => {A,B}
                        if (child.right_before_next()) {
                            child.activate_state(node_state::last_branching_after_left_end);
                        } else {
                            child.activate_state(node_state::branching_after_left_end);
                        }
                        break;
                } case node_state::non_branching_left_only: { // D => {G}
                        child.activate_state(node_state::non_branching_after_left);
                        break;
                } case node_state::non_branching_after_left: {[[fallthrough]];
                } case node_state::non_branching_including_left: { // {G,H} => {A,B} | {G}
                        if (child.get_right() != child.sink() && libjst::position(*child.get_right()).is_left_end()) {
                            if (child.right_before_next()) { // => {B}
                                child.activate_state(node_state::last_branching_after_left_end);
                            } else { // => {A}
                                child.activate_state(node_state::branching_after_left_end);
                            }
                        } else { // => {G}
                            child.activate_state(node_state::non_branching_after_left);
                        }
                        break;
                } case node_state::branching_after_left_end: { // A => {A/B}
                            if (child.right_before_next()) { // => {B}
                                child.activate_state(node_state::last_branching_after_left_end);
                            } else { // => {A}
                                child.activate_state(node_state::branching_after_left_end);
                            }
                            break;
                } case node_state::branching_after_left_begin: { // E => {E/F}
                            if (child.right_before_next()) { // => {F}
                                child.activate_state(node_state::last_branching_after_left_begin);
                            } else { // => {E}
                                child.activate_state(node_state::branching_after_left_begin);
                            }
                            break;
                }
                default: /*no-op*/;
            }
            return child;
        }

        constexpr bool right_before_next() const noexcept {
            return std::ranges::distance(get_right(), get_next()) == 1;
        }

        static constexpr rcs_store_node_base & as_base(derived_t & derived) noexcept {
            return static_cast<rcs_store_node_base &>(derived);
        }

        static constexpr rcs_store_node_base const & as_base(derived_t const & derived) noexcept {
            return static_cast<rcs_store_node_base const &>(derived);
        }

        static constexpr derived_t & as_derived(rcs_store_node_base & base) noexcept {
            return static_cast<derived_t &>(base);
        }

        static constexpr derived_t const & as_derived(rcs_store_node_base const & base) noexcept {
            return static_cast<derived_t const &>(base);
        }

        constexpr friend bool operator==(rcs_store_node_base const & lhs, rcs_store_node_base const & rhs) noexcept {
            [[maybe_unused]] bool s1 = static_cast<node_descriptor const &>(lhs) == static_cast<node_descriptor const &>(rhs);
            [[maybe_unused]] bool s2 = lhs._left_variant == rhs._left_variant;
            [[maybe_unused]] bool s3 = lhs._right_variant == rhs._right_variant;
            return  s1 && s2 && s3;
                //    std::tie(lhs._left_variant, lhs._right_variant) == std::tie(rhs._left_variant, rhs._right_variant);
        }
    };
}  // namespace libjst
