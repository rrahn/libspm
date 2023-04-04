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

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>
#include <libjst/sequence_tree/rcs_node_traits.hpp>
#include <libjst/variant/breakpoint.hpp>

namespace libjst
{
    template <typename base_tree_t>
    class partial_tree {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using variant_type = typename rcs_node_traits<base_node_type>::variant_type;
        using breakpoint_value_type = typename breakpoint::value_type;

        class node_impl;
        class sink_impl;

        base_tree_t _wrappee{};
        variant_type _left_bound{};
        variant_type _right_bound{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr partial_tree() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<wrapped_tree_t, partial_tree> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr partial_tree(wrapped_tree_t && wrappee,
                                        breakpoint_value_type root_position,
                                        size_t count) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee}
        {
            _left_bound = *std::ranges::begin(_wrappee.data().variants());
            _right_bound = _left_bound;

            breakpoint_value_type sink_position =
                std::min<size_t>(std::ranges::size(_wrappee.data().source()), root_position + count);
            libjst::position(_left_bound) = breakpoint{root_position, breakpoint_end::right};
            libjst::position(_right_bound) = breakpoint{sink_position, breakpoint_end::left};
        }
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee), _left_bound, _right_bound};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(_wrappee)};
        }

        constexpr base_tree_t const & base() const noexcept {
            return _wrappee;
        }
   };

    template <typename wrapped_tree_t, std::integral offset_t, std::integral count_t>
    partial_tree(wrapped_tree_t &&, offset_t, count_t) -> partial_tree<wrapped_tree_t>;

    template <typename base_tree_t>
    class partial_tree<base_tree_t>::node_impl : public base_node_type {
    private:

        using variant_iterator = typename rcs_node_traits<base_node_type>::variant_iterator;
        using variant_reference = std::iter_reference_t<variant_iterator>;

        friend partial_tree;

        enum state {
            left_bound,
            right_bound,
            regular,
        };

        variant_type const * _left_bound{};
        variant_type const * _right_bound{};
        state _left_state{state::left_bound};
        state _right_state{state::regular};

        explicit constexpr node_impl(base_node_type && base_node,
                                     variant_type const & left_bound,
                                     variant_type const & right_bound) noexcept :
            base_node_type{std::move(base_node)},
            _left_bound{std::addressof(left_bound)},
            _right_bound{std::addressof(right_bound)}
        {
            // now set the variants to the correct positions! -> would it be possible to add a breakpoint layer that does this for us?
            // we do not have a root node, but can emulate one!
            // we are somewhere between left node and right node.
            // left variant => largest variant before root_breakpoint
            auto const & variants = base_node_type::rcs_store().variants();
            auto initial_it = std::ranges::lower_bound(std::ranges::next(variants.begin()),
                                                       variants.end(),
                                                       libjst::position(*_left_bound),
                                                       std::ranges::less{},
                                                       [] (auto const & var) { return libjst::position(var); });
            assert(initial_it != variants.begin());
            base_node_type::set_left(std::ranges::prev(initial_it));
            base_node_type::set_right(initial_it);
            base_node_type::set_next(base_node_type::next_variant_after(initial_it));
            // initialise left state
            if (libjst::position(*_left_bound).value() == 0)
                _left_state = state::regular;
            else
                _left_state = state::left_bound;

            // initialise right state
            if (libjst::position(base_node_type::right_variant()) < libjst::position(*_right_bound))
                _right_state = state::regular;
            else
                _right_state = state::right_bound;
        }

        explicit constexpr node_impl(base_node_type && base_node,
                                     variant_type const * left_bound,
                                     variant_type const * right_bound,
                                     state left_state,
                                     state right_state) noexcept :
            base_node_type{std::move(base_node)},
            _left_bound{left_bound},
            _right_bound{right_bound},
            _left_state{left_state},
            _right_state{right_state}
        {}

    public:

        constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit<false>(base_node_type::next_ref());
        }

    protected:

        constexpr breakpoint left_breakpoint() const noexcept {
            return select_breakpoint(_left_state, base_node_type::left_breakpoint());
        }

        constexpr breakpoint right_breakpoint() const noexcept {
            return select_breakpoint(_right_state, base_node_type::right_breakpoint());
        }

        constexpr variant_reference left_variant() const noexcept {
            switch (_left_state) {
                case state::left_bound: return *_left_bound;
                case state::right_bound: return *_right_bound;
                default: return base_node_type::left_variant();
            }
        }

        constexpr variant_reference right_variant() const noexcept {
            switch (_right_state) {
                case state::left_bound: return *_left_bound;
                case state::right_bound: return *_right_bound;
                default: return base_node_type::right_variant();
            }
        }

        constexpr bool is_branching() const noexcept {
            return !reached_right_bound() && base_node_type::is_branching();
        }

        constexpr bool is_nil() const noexcept {
            return base_node_type::is_nil() || reached_right_bound();
        }
    private:

        constexpr bool reached_right_bound() const noexcept {
            bool not_alternate = !base_node_type::on_alternate_path();
            bool over_right_bound = libjst::position(*_right_bound) <= base_node_type::right_breakpoint();
            return not_alternate && over_right_bound;
        }

        constexpr breakpoint select_breakpoint(state bound_state, breakpoint alt_breakpoint) const noexcept {
            switch (bound_state) {
                case state::left_bound: return breakpoint{(uint32_t)libjst::position(*_left_bound), breakpoint_end::left};
                case state::right_bound: return breakpoint{(uint32_t)libjst::position(*_right_bound), breakpoint_end::left};
                default: return alt_breakpoint;
            }
        }

        template <bool is_alt>
        constexpr std::optional<node_impl> visit(auto maybe_child) const {
            if (maybe_child) {
                state child_right_state = (is_alt) ? state::regular : _right_state;
                return node_impl{std::move(*maybe_child), _left_bound, _right_bound, state::regular, child_right_state};
            } else {
                return std::nullopt;
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return lhs.reached_right_bound() || static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class partial_tree<base_tree_t>::sink_impl {
    private:
        friend partial_tree;

        using base_sink_type = libjst::tree_sink_t<base_tree_t>;
        base_sink_type _base_sink{};

        constexpr explicit sink_impl(base_sink_type base_sink) : _base_sink{std::move(base_sink)}
        {}

        friend bool operator==(sink_impl const & lhs, base_node_type const & rhs) noexcept {
            return lhs._base_sink == rhs;
        }

    public:
        sink_impl() = default;
    };
}  // namespace libjst
