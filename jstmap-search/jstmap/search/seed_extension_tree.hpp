// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides extension base tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/variant/concept.hpp>

#include <jstmap/global/match_position.hpp>

namespace jstmap
{
    template <typename wrapped_tree_t>
    class seed_extension_tree
    {
    private:
        using base_node_type = libjst::tree_node_t<wrapped_tree_t>;
        using sink_type = libjst::tree_sink_t<wrapped_tree_t>;
        using base_cargo_type = libjst::tree_label_t<wrapped_tree_t>;
        using difference_type = std::ptrdiff_t;

        // generate local position!

        class node_impl;
        class cargo_impl;

        wrapped_tree_t _wrappee{};
        base_node_type _base{};
        difference_type _low_position{};
        difference_type _max_label_size{};


    public:

        template <typename wrappee_t, std::unsigned_integral extension_size_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, seed_extension_tree>) &&
                      std::constructible_from<wrapped_tree_t, wrappee_t>
        seed_extension_tree(wrappee_t && wrappee, match_position start, extension_size_t max_label_size) noexcept :
            _wrappee{(wrappee_t &&) wrappee}
        {
        // case a) start on reference path
            // - seek to the current position in the root node // that is we know the variant index
            auto tmp_tree = _wrappee | libjst::seek();
            _base = tmp_tree.seek(start.tree_position).base(); // from seeking to the node!
            auto base_cargo = *_base;
            // - compute the label distance until next breakend.
            assert(std::ranges::ssize(base_cargo.path_sequence()) >= start.label_offset);
            difference_type _distance_to_high = std::ranges::ssize(base_cargo.path_sequence()) - start.label_offset;
            _low_position = libjst::position(_base.high_boundary()) - _distance_to_high;
            assert(static_cast<difference_type>(libjst::position(_base.low_boundary())) <= _low_position);
            // - update remaining max label size of root node.
            _max_label_size = static_cast<difference_type>(max_label_size) - _distance_to_high;
        }

        constexpr node_impl root() const noexcept {
            return node_impl{_base, _low_position, _max_label_size};
        }

        constexpr sink_type sink() const noexcept {
            return _wrappee.sink();
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
    };

    template <typename wrapped_tree_t>
    class seed_extension_tree<wrapped_tree_t>::node_impl : public base_node_type
    {
    private:
        friend seed_extension_tree;

        // using base_low_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().low_boundary())>;
        // using base_high_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().high_boundary())>;

        friend seed_extension_tree;

        difference_type _low_position{};
        difference_type _remaining_label_size{};

        explicit constexpr node_impl(base_node_type base_node,
                                     difference_type low_position,
                                     difference_type remaining_label_size) noexcept :
            base_node_type{std::move(base_node)},
            _low_position{low_position},
            _remaining_label_size{remaining_label_size}
        {}

    public:

        // TODO: Overload the low position type!

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (is_leaf())
                return std::nullopt;
            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (is_leaf())
                return std::nullopt;
            return visit<false>(base_node_type::next_ref());
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

    private:

        constexpr bool is_leaf() const noexcept {
            return _remaining_label_size <= 0;
        }

        template <bool is_alt_node, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                // root node already subtracted the offset.
                // was called on alt node and parent not on alternate path
                difference_type next_span = libjst::position(maybe_child->high_boundary()) -
                                            libjst::position(maybe_child->low_boundary());
                if constexpr (is_alt_node) {
                    auto delta = *(maybe_child->low_boundary());
                    next_span += libjst::effective_size(delta) - std::ranges::ssize(libjst::alt_sequence(delta));
                }
                return node_impl{std::move(*maybe_child), _low_position, _remaining_label_size - next_span};
            } else {
                return std::nullopt;
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename wrapped_tree_t>
    class seed_extension_tree<wrapped_tree_t>::cargo_impl : public base_cargo_type
    {
    private:
        friend seed_extension_tree;

        node_impl const * _node{};
        constexpr explicit cargo_impl(node_impl const * node) noexcept :
            base_cargo_type{*static_cast<base_node_type const &>(*node)},
            _node{node}
        {
        }

    public:

        constexpr auto sequence() const noexcept {
            assert(_node != nullptr);

            return base_cargo_type::sequence(get_low_position(), get_high_position());
        }

        constexpr auto path_sequence() const noexcept {
            assert(_node != nullptr);
            return base_cargo_type::sequence(0, get_high_position());
        }

    protected:
        using base_cargo_type::sequence;

    private:

        constexpr difference_type get_low_position() const noexcept {
            difference_type low_position = libjst::position(_node->low_boundary());
            return std::max(low_position, _node->_low_position);
        }

        constexpr difference_type get_high_position() const noexcept {
            difference_type high_position = libjst::position(_node->high_boundary());
            return std::min(high_position, high_position + _node->_remaining_label_size);
        }
    };

    namespace _tree_adaptor {
        struct _extend_from
        {
            template <typename tree_t, std::unsigned_integral extension_size_t>
            constexpr auto operator()(tree_t && tree, match_position start_position, extension_size_t extension_size) const
                noexcept(std::is_nothrow_constructible_v<seed_extension_tree<std::remove_reference_t<tree_t>>>)
                -> seed_extension_tree<std::remove_reference_t<tree_t>>
            {
                using adapted_tree_t = seed_extension_tree<std::remove_reference_t<tree_t>>;
                return adapted_tree_t{(tree_t &&)tree, std::move(start_position), extension_size};
            }

            template <std::unsigned_integral extension_size_t>
            constexpr auto operator()(match_position start_position, extension_size_t extension_size) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _extend_from, match_position, extension_size_t>)
                -> jst::contrib::closure_result_t<_extend_from, match_position, extension_size_t>
            {
                return jst::contrib::make_closure(_extend_from{}, std::move(start_position), extension_size);
            }
        };
        inline constexpr _extend_from extend_from{};
    } // namespace _tree_adaptor

    using _tree_adaptor::extend_from;
}  // namespace jstmap
