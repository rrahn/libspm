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

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/node_descriptor.hpp>
#include <libjst/variant/breakpoint.hpp>

namespace libjst
{
    template <typename base_tree_t>
    class merge_tree_impl {
    private:
        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;
        using base_cargo_type = libjst::tree_label_t<base_tree_t>;

        class node_impl;
        class cargo_impl;

        base_tree_t _wrappee{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr merge_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
            requires (!std::same_as<wrapped_tree_t, merge_tree_impl> &&
                      std::constructible_from<base_tree_t, wrapped_tree_t>)
        explicit constexpr merge_tree_impl(wrapped_tree_t && wrappee) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            base_node_type base_root = libjst::root(_wrappee);
            auto root_low = base_root.low_boundary();
            return node_impl{std::move(base_root), std::move(root_low)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }

        constexpr auto const & data() const noexcept {
            return _wrappee.data();
        }
   };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::node_impl : public base_node_type {
    public:

        using low_position_type = std::remove_cvref_t<decltype(std::declval<base_node_type const &>().low_boundary())>;

    private:

        friend merge_tree_impl;

        low_position_type _low_boundary{};

        explicit constexpr node_impl(base_node_type && base_node, low_position_type cached_low) noexcept :
            base_node_type{std::move(base_node)},
            _low_boundary{std::move(cached_low)}
        {}

    public:

        constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            return visit_next<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            return visit_next<false>(base_node_type::next_ref());
        }

        constexpr low_position_type const & low_boundary() const {
            return _low_boundary;
        }

        constexpr cargo_impl operator*() const noexcept {
            return cargo_impl{this};
        }

    protected:

        template <typename breakend_site_t>
        constexpr void reset_low(breakend_site_t new_low) {
            base_node_type::reset_low(std::move(new_low));
            _low_boundary = base_node_type::low_boundary();
            extend();
        }

    private:

        template <bool is_alt_child>
        constexpr std::optional<node_impl> visit_next(auto maybe_child) const {
            if (maybe_child) {
                low_position_type cached_low = maybe_child->low_boundary();
                node_impl new_child{std::move(*maybe_child), std::move(cached_low)};
                new_child.extend();
                return new_child;
            } else {
                return std::nullopt;
            }
        }

        constexpr void extend() {
            while (!base_node_type::high_boundary().is_low_end()) {
                if (auto successor = base_node_type::next_ref(); successor) {
                    static_cast<base_node_type &>(*this) = std::move(*successor);

                } else {
                    break;
                }
            }
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept
        {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class merge_tree_impl<base_tree_t>::cargo_impl : public base_cargo_type {
    private:
        using position_type = typename base_node_type::position_type;

        friend merge_tree_impl;

        node_impl const * _node{};

        explicit constexpr cargo_impl(node_impl const * node) noexcept :
            base_cargo_type{*static_cast<base_node_type const &>(*node)},
            _node{node}
        {}

    public:

        cargo_impl() = default;

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
        inline constexpr struct _merge
        {
            template <typename covered_tree_t, typename ...args_t>
            constexpr auto operator()(covered_tree_t && tree, args_t &&... args) const
                noexcept(std::is_nothrow_constructible_v<
                            merge_tree_impl<std::remove_reference_t<covered_tree_t>>, args_t...>)
                -> merge_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>
            {
                using adapted_tree_t = merge_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>;
                return adapted_tree_t{(covered_tree_t &&)tree, (args_t &&)args...};
            }

            template <typename ...args_t>
            constexpr auto operator()(args_t &&... args) const
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<jst::contrib::make_closure>, args_t...>)
                -> jst::contrib::closure_result_t<_merge, args_t...>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_merge{}, (args_t &&)args...);
            }
        } merge{};
    } // namespace _tree_adaptor

    using _tree_adaptor::merge;
}  // namespace libjst
