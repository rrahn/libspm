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

#include <libjst/utility/closure_object.hpp>

#include <libjst/sequence_tree/breakend_site_trimmed.hpp>
#include <libjst/sequence_tree/concept.hpp>
namespace libjst
{
    template <typename base_tree_t>
    class trim_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using base_cargo_type = libjst::tree_label_t<base_tree_t>;
        using difference_type = std::ptrdiff_t;

        class node_impl;
        class cargo_impl;

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

        using base_low_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().low_boundary())>;
        using base_high_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().high_boundary())>;

        friend trim_tree_impl;

        difference_type _max_branch_size{};

        explicit constexpr node_impl(base_node_type base_node, difference_type max_branch_size) noexcept :
            base_node_type{std::move(base_node)},
            _max_branch_size{max_branch_size}
        {}

    public:

        using low_position_type = base_low_position_type;
        using high_position_type = breakend_site_trimmed<base_high_position_type>;

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

        constexpr high_position_type high_boundary() const {
            base_high_position_type base_high = base_node_type::high_boundary();
            if (this->on_alternate_path()) {
                using position_value_t = typename high_position_type::position_value_type;
                // TODO: Assumes not left extended! or we need to add the left extension below!
                position_value_t high_position = libjst::position(base_high);
                assert(static_cast<difference_type>(high_position) + _max_branch_size > 0);
                return high_position_type{std::move(base_high),
                                          static_cast<position_value_t>(high_position + _max_branch_size)};
            } else {
                return high_position_type{std::move(base_high)};
            }
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

    private:

        constexpr bool is_leaf() const noexcept {
            return _max_branch_size <= 0;
        }

        template <bool is_alt_node, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                if (is_alt_node && !this->on_alternate_path()) {
                    return make_alternate_subtree(std::move(*maybe_child));
                } else if (this->on_alternate_path()) {
                    return branch_off_further<is_alt_node>(std::move(*maybe_child));
                } else { // nothing to spawn - remain in the reference branch.
                    return node_impl{std::move(*maybe_child), _max_branch_size};
                }
            } else {
                return std::nullopt;
            }
        }

        constexpr node_impl make_alternate_subtree(base_node_type && base_child) const noexcept {
            auto delta = *(base_child.low_boundary());
            difference_type diff =
                libjst::position(base_child.high_boundary()) - libjst::position(base_child.low_boundary()) +
                libjst::effective_size(delta) - std::ranges::ssize(libjst::alt_sequence(delta));
            return node_impl{std::move(base_child), _max_branch_size - diff};
        }

        template <bool is_alt_node>
        constexpr node_impl branch_off_further(base_node_type && base_child) const noexcept {
            difference_type child_remaining = _max_branch_size;
            if constexpr (is_alt_node) {
                child_remaining -= std::ranges::ssize(libjst::alt_sequence(*base_child.low_boundary()));
            } else {
                child_remaining -= (libjst::position(base_child.high_boundary()) -
                                    libjst::position(base_child.low_boundary()));
            }

            return node_impl{std::move(base_child), child_remaining};
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class trim_tree_impl<base_tree_t>::cargo_impl : public base_cargo_type {
    private:
        friend trim_tree_impl;

        node_impl const * _node{};

        explicit constexpr cargo_impl(node_impl const * node) noexcept :
            base_cargo_type{*static_cast<base_node_type const &>(*node)},
            _node{node}
        {}

    public:
        constexpr cargo_impl() = default;

        constexpr auto sequence() const noexcept {
            assert(_node != nullptr);
            return base_cargo_type::sequence(libjst::position(_node->low_boundary()),
                                             libjst::position(_node->high_boundary()));
        }

        constexpr auto path_sequence() const noexcept {
            assert(_node != nullptr);
            return base_cargo_type::sequence(0, libjst::position(_node->high_boundary()));
        }

    protected:
        using base_cargo_type::sequence;
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
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<libjst::make_closure>, _trim, branch_size_t>)
                -> libjst::closure_result_t<_trim, branch_size_t>
            {
                return libjst::make_closure(_trim{}, branch_size);
            }
        };
        inline constexpr _trim trim{};
    } // namespace _tree_adaptor

    using _tree_adaptor::trim;
}  // namespace libjst
