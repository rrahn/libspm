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

    private:
        rcs_store_type const * _rcs_store{};
        variant_iterator _left_variant{};
        variant_iterator _right_variant{};
        variant_iterator _next_variant{};

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

            node_descriptor::set_reference();
            // initial state?
            node_descriptor::set_first_breakpoint_id(node_descriptor_id::first_right);
            if (_right_variant != sink()) {
                node_descriptor::set_second_breakpoint_id(node_descriptor_id::second_left);
            } else {
                node_descriptor::set_second_breakpoint_id(node_descriptor_id::second_right);
            }

            set_next(next_variant_after(_right_variant));
        }

        // only from branching node!
        constexpr optional_node_type visit_next_alt() const noexcept {
            if (is_ref_node() && is_branching()) {
                derived_t child{as_derived(*this)};
                child.set_left(get_right());
                child.set_alternate();
                child.set_first_breakpoint_id(node_descriptor_id::first_left);
                child.set_second_breakpoint_id(node_descriptor_id::second_first_right);
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
                child.set_reference();
                child.set_first_breakpoint_id(node_descriptor_id::first_right);
                if (child.get_right() != sink() && position(*child.get_right()).is_left_end()) {
                    child.set_second_breakpoint_id(node_descriptor_id::second_left);
                } else {
                    child.set_second_breakpoint_id(node_descriptor_id::second_right);
                }
            }
            return child;
        }

        constexpr bool is_ref_node() const noexcept {
            return node_descriptor::from_reference();
        }

        constexpr bool is_alt_node() const noexcept {
            return node_descriptor::from_alternate();
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
        constexpr breakpoint_type left_breakpoint() const {
            if (get_left() == sink())
                return std::ranges::size(rcs_store().source());

            auto bounded_right_position = [&] (variant_iterator const it) {
                if (it == sink()) {
                    return breakpoint{static_cast<uint32_t>(std::ranges::size(rcs_store().source()))};
                } else {
                    return libjst::position(*it);
                }
            };

            // return (is_variant_link()) ? libjst::right_breakpoint(left_var) : libjst::left_breakpoint(left_var);
            using seqan3::operator&;
            return ((node_descriptor::get_first_breakpoint_id() & node_descriptor_id::first_left) == node_descriptor_id::first_left)
                        ? libjst::left_breakpoint(*get_left())
                        : std::min(libjst::right_breakpoint(*get_left()), bounded_right_position(get_right()));
        }

        /*!\brief Returns the breakpoint of the right end of the current node.
         *
         * The correct position of the current node depends on the node type. If `this` is not a branching node, i.e.
         * libjst::rcs_store_node_base::is_branching_node() evaluates to `false`, this member returns the position of
         * the right variant, otherwise it returns the end position of the left variant.
         *
         * \returns Source position of the right end of the current node.
         */
        constexpr breakpoint_type right_breakpoint() const {

            auto bounded_right_breakpoint = [&] (variant_iterator const it) {
                if (it == sink()) {
                    return breakpoint{static_cast<uint32_t>(std::ranges::size(rcs_store().source()))};
                } else {
                    return libjst::right_breakpoint(*it);
                }
            };

            auto bounded_left_breakpoint = [&] (variant_iterator const it) {
                if (it == sink()) {
                    return breakpoint{static_cast<uint32_t>(std::ranges::size(rcs_store().source()))};
                } else {
                    return libjst::left_breakpoint(*it);
                }
            };

            using seqan3::operator&;
            if ((node_descriptor::get_second_breakpoint_id() & node_descriptor_id::second_first_right) != node_descriptor_id::nil) {
                return bounded_right_breakpoint(get_left());
            } else if ((node_descriptor::get_second_breakpoint_id() & node_descriptor_id::second_right) != node_descriptor_id::nil) {
                return bounded_right_breakpoint(get_right());
            } else {
                return bounded_left_breakpoint(get_right());
            }
            // return (node_descriptor::get_second_breakpoint_id() == node_descriptor_id::first_right)
            //             ? libjst::left_breakpoint(*left_variant())
            //             : libjst::right_breakpoint(*left_variant());
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
            return get_right() != sink() && node_descriptor::is_branching();
        }

        constexpr variant_iterator sink() const noexcept {
            return std::ranges::end(rcs_store().variants());
        }

    private:

        constexpr bool is_leaf_of_alternate_subtree() const noexcept {
            return node_descriptor::on_alternate_path() &&
                   get_right() == sink() &&
                   node_descriptor::get_second_breakpoint_id() != node_descriptor_id::second_first_right;
        }

        constexpr variant_iterator get_next() const noexcept {
            return _next_variant;
        }
        // constexpr bool is_variant_link() const noexcept {
        //     return libjst::left_breakpoint(*left_variant()) != libjst::left_breakpoint(*right_variant());
        // }

        // constexpr void toggle_branching_state() noexcept {
        //     assert(node_descriptor::from_reference());
        //     if (right_variant() != sink() && ) {
        //         node_descriptor::set_branching();
        //     } else {
        //         node_descriptor::unset_branching();
        //     }
        // }

        constexpr variant_iterator find_next_valid_right_variant() const noexcept {
            variant_iterator next_right = get_right();
            // if (is_ref_node()) {
            //     next_right = std::ranges::next(next_right, static_cast<difference_type>(!node_descriptor::is_branching()), sink());
            // } else {
                assert(is_alt_node());
                auto const min_ref_position = libjst::right_breakpoint(*next_right);
                next_right = std::ranges::find_if_not(std::ranges::next(next_right), sink(), [&] (auto && var) {
                    return libjst::left_breakpoint(var) < min_ref_position;
                });
                assert(next_right != get_right());
            // }
            return next_right;
        }

        static constexpr derived_t visit_next_ref_impl(derived_t child) noexcept {
            auto parent_node_state = get_ref_node_state(child);
            // Move through the states.
            switch(parent_node_state) {
                case node_state::A: [[fallthrough]];
                case node_state::E:
                    child.set_left(child.get_right());
                    child.set_right(std::ranges::next(child.get_right()));
                    break;
                case node_state::B: [[fallthrough]];
                case node_state::F: [[fallthrough]];
                case node_state::G: [[fallthrough]];
                case node_state::H:
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
                  case node_state::B: { [[fallthrough]];
                } case node_state::F: {
                    child.set_first_breakpoint_id(node_descriptor_id::first_left);
                    if (libjst::right_breakpoint(*child.get_left()) < bounded_left_breakpoint(child.get_right())) {
                        child.set_second_breakpoint_id(node_descriptor_id::second_first_right); // C/D
                    } else if (bounded_breakpoint(child.get_right()).is_left_end()) {
                        child.set_second_breakpoint_id(node_descriptor_id::second_left); // E/F
                    } else {
                        child.set_second_breakpoint_id(node_descriptor_id::second_right); // H
                    }
                    break;
                } case node_state::C: { // => {A,B}
                        child.set_first_breakpoint_id(node_descriptor_id::first_right);
                        child.set_second_breakpoint_id(node_descriptor_id::second_left);
                        break;
                } case node_state::D: { // => {G}
                        child.set_first_breakpoint_id(node_descriptor_id::first_right);
                        child.set_second_breakpoint_id(node_descriptor_id::second_right);
                        break;
                } case node_state::G: {[[fallthrough]];
                } case node_state::H: { // => {A,B} | {G}
                        child.set_first_breakpoint_id(node_descriptor_id::first_right);
                        if (child.get_right() != child.sink() && libjst::position(*child.get_right()).is_left_end()) {
                            child.set_second_breakpoint_id(node_descriptor_id::second_left); // A,B
                        } else {
                            child.set_second_breakpoint_id(node_descriptor_id::second_right); // G
                        }
                        break;
                } //case node_state::A: { [[fallthrough]]; // either in A or in B -> same properties
                //} case node_state::E: { return child; // either in E or in F -> same properties
                //}
                default: /*no-op*/;
            }
            return child;
        }

        static constexpr node_state get_ref_node_state(derived_t const & node) noexcept {
            using seqan3::operator&;
            using seqan3::operator|;
            assert((node.get_first_breakpoint_id() & (node_descriptor_id::first_left | node_descriptor_id::first_right)) !=
                   node_descriptor_id::nil);

            if ((node.get_first_breakpoint_id() & node_descriptor_id::first_left) != node_descriptor_id::nil) {
                if ((node.get_second_breakpoint_id() & node_descriptor_id::second_first_right) != node_descriptor_id::nil) {
                    [[maybe_unused]] bool tmpB = node.get_right() == node.sink();
                    [[maybe_unused]] bool tmpA = node.get_right() != node.sink();
                    return (tmpA && libjst::position(*node.get_right()).is_left_end())
                                ? node_state::C
                                : node_state::D;
                } else { // node.get_second_breakpoint_id() & node_descriptor_id::second_right -> true
                    assert((node.get_second_breakpoint_id() & (node_descriptor_id::second_left | node_descriptor_id::second_right)) != node_descriptor_id::nil);
                    if ((node.get_second_breakpoint_id() & node_descriptor_id::second_left) != node_descriptor_id::nil) {
                        return (std::ranges::distance(node.get_right(), node.get_next()) > 1)
                                ? node_state::E
                                : node_state::F;
                    } else {
                        return node_state::H;
                    }
                }
            } else { // node.get_first_breakpoint_id() & node_descriptor_id::first_right -> true
                assert((node.get_first_breakpoint_id() & node_descriptor_id::first_right) != node_descriptor_id::nil);
                if ((node.get_second_breakpoint_id() & node_descriptor_id::second_left) != node_descriptor_id::nil) {
                    return (std::ranges::distance(node.get_right(), node.get_next()) > 1)
                                ? node_state::A
                                : node_state::B;
                } else {
                    return node_state::G;
                }
            }
        }

        // States of parent node when calling next_ref()
        // * node_descriptor::is_reference()
        //      * node_descriptor::is_branching()
        //          * A: (s-j) > 1, j = left end  => create branching ref child: ++j                         | i < j < s, j = left end -> {A,B}
        //          * B: (s-j) == 1, j = left end => create non-branching ref child: i = j, j = s, s get nxt | i < j < s               -> {C, D}
        //      * !node_descriptor::is_branching()
        //          * i = left end
        //              * C: i < j < s, j = left end  => create branching ref child:                               | i < j < s, j = left end -> {A, B}
        //              * D: i < j < s, j = right end => create non-branching ref child: i = j, j = s, s get next  | i < j < s, i = right end -> {E, F}
        //          * i = right end
        //              * E: i < j < s, j = left end  => create branching ref child:                              | i < j < s, j = left end  -> {A, B}
        //              * F: i < j < s, j = right end => create non-branching ref child: i = j, j = s, s get next | i < j < s, i = right end -> {E, F}
        // * node_descriptor::is_alternate()
        //   * G: i = j < s, i = j = left, => create ref child: i = j, j = first_after_right_end_of(i), s = get nxt(j) | i < j < s -> {A, B, E, F}

        // label states:
        // A -> lbl[pr(i), pl(j)]
        // B -> lbl[pr(i), pl(j)]
        // C -> lbl[pl(i), min(pr(i), pl(j))]
        // D -> lbl[pl(i), min(pr(i), pr(j))]
        // E -> lbl[pr(i), pl(j)]
        // F -> lbl[pr(i), pr(j)]
        // G -> lbl[pl(i), pr(j)]

        // position states:
        // A -> pr(i) <= pl(j)

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
            return static_cast<node_descriptor const &>(lhs) == static_cast<node_descriptor const &>(rhs) &&
                   std::tie(lhs._left_variant, lhs._right_variant) == std::tie(rhs._left_variant, rhs._right_variant);
        }
    };
}  // namespace libjst
