// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides trim tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/breakend_site_trimmed.hpp>
#include <libjst/sequence_tree/concept.hpp>
namespace libjst
{
    template <typename base_tree_t>
    class trim_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using difference_type = std::ptrdiff_t;

        class node_impl;

        base_tree_t _wrappee{};
        difference_type _max_branch_size{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr trim_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<std::remove_cvref_t<wrapped_tree_t>, trim_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr trim_tree_impl(wrapped_tree_t && wrappee, std::size_t const max_branch_size) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee},
            _max_branch_size{static_cast<difference_type>(max_branch_size)}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee), _max_branch_size};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
   };

    template <typename base_tree_t>
    class trim_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:
        using base_position_type = typename base_node_type::position_type;

        friend trim_tree_impl;

        difference_type _max_branch_size{};
        difference_type _remaining_branch_size{};

        explicit constexpr node_impl(base_node_type base_node, difference_type max_branch_size) noexcept :
            node_impl{std::move(base_node), max_branch_size, max_branch_size}
        {}

        explicit constexpr node_impl(base_node_type && base_node,
                                     difference_type max_branch_size,
                                     difference_type remaining_branch_size) noexcept :
            base_node_type{std::move(base_node)},
            _max_branch_size{max_branch_size},
            _remaining_branch_size{remaining_branch_size}
        {}

    public:

        using position_type = breakend_site_trimmed<base_position_type>;

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

        constexpr libjst::node_label_t<base_node_type> operator*() const noexcept {
            return *static_cast<base_node_type const &>(*this);
        }

        // If we cut the tree based on the length we need to specify different end position.
        // So far we have used it through the values.
        constexpr position_type low_boundary() const {
            return position_type{base_node_type::low_boundary()};
        }

        constexpr position_type high_boundary() const {
            using position_value_t = typename position_type::position_value_type;
            base_position_type const & base_boundary = base_node_type::high_boundary();
            if (this->on_alternate_path()) {
                assert(static_cast<difference_type>(libjst::position(base_boundary)) + _remaining_branch_size > 0);
                return position_type{base_boundary, static_cast<position_value_t>(libjst::position(base_boundary) +
                                                                                  _remaining_branch_size)};
                        // std::min<difference_type>(static_cast<difference_type>(original_bp),
                        //                           static_cast<difference_type>(original_bp) + _remaining_branch_size));
            } else {
                return position_type{base_boundary};
            }
        }
    private:

        constexpr bool is_leaf() const noexcept {
            return _remaining_branch_size <= 0;
        }

        template <bool is_alt_node, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                if (is_alt_node && !this->on_alternate_path()) { // parent not on alternate path but child is
                    return branch_off_new(std::move(*maybe_child));
                } else if (base_node_type::on_alternate_path()) { // parent on alternate path
                    return branch_off_further<is_alt_node>(std::move(*maybe_child));
                } else { // nothing to spawn - remain in the reference branch.
                    return node_impl{std::move(*maybe_child), _max_branch_size, _remaining_branch_size};
                }
            } else {
                return std::nullopt;
            }
        }

        // branches off the reference path
        constexpr node_impl branch_off_new(base_node_type && base_child) const noexcept {
            difference_type child_branch_size = _max_branch_size +
                                std::ranges::ssize(libjst::alt_sequence(*(base_node_type::high_boundary())));
            return node_impl{std::move(base_child), child_branch_size, _remaining_branch_size};
        }

        template <bool is_alt_node>
        constexpr node_impl branch_off_further(base_node_type && base_child) const noexcept {
            difference_type child_remaining = _remaining_branch_size;
            if constexpr (is_alt_node) {
                child_remaining -= std::ranges::ssize(libjst::alt_sequence(*base_child.low_boundary()));
            } else {
                child_remaining -= (libjst::position(base_child.high_boundary()) -
                                    libjst::position(base_child.low_boundary()));
            }

            return node_impl{std::move(base_child), _max_branch_size, child_remaining};
            // std::cout << "Further Span: " << new_child.breakend_span<is_alt_node>() << "\n";
            // if (new_child.breakend_span<is_alt_node>() > 1000) {
            //     std::cout << "Stop\n";
            // }
            // new_child._remaining_branch_size -= new_child.breakend_span<is_alt_node>();
            // return new_child;
        }

        // template <bool is_alt_node>
        // constexpr difference_type breakend_span() const noexcept {
        //     if constexpr (is_alt_node) {
        //         return std::ranges::ssize(libjst::alt_sequence(base_node_type::left_variant())); //TODO
        //     } else {
        //         return base_node_type::high_breakend() - base_node_type::low_breakend();
        //     }
        // }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    namespace _tree_adaptor {
        struct _trim
        {
            template <typename tree_t, std::unsigned_integral branch_size_t>
            constexpr auto operator()(tree_t && tree, branch_size_t const branch_size) const
                noexcept(std::is_nothrow_constructible_v<trim_tree_impl<std::remove_reference_t<tree_t>>>)
                -> trim_tree_impl<std::remove_reference_t<tree_t>>
            {
                using adapted_tree_t = trim_tree_impl<std::remove_reference_t<tree_t>>;
                return adapted_tree_t{(tree_t &&)tree, branch_size};
            }

            template <std::unsigned_integral branch_size_t>
            constexpr auto operator()(branch_size_t const branch_size) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, _trim, branch_size_t>)
                -> jst::contrib::closure_result_t<_trim, branch_size_t>
            {
                return jst::contrib::make_closure(_trim{}, branch_size);
            }
        };
        inline constexpr _trim trim{};
    } // namespace _tree_adaptor

    using _tree_adaptor::trim;
}  // namespace libjst
