// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides k_depth tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/copyable_box.hpp>

// #include <libjst/sequence_tree/rcs_node_traits.hpp>
#include <libjst/sequence_tree/concept.hpp>

namespace libjst
{
    template <typename base_tree_t>
        // requires covered tree
    class k_depth_tree_impl {
    private:
        using wrappee_t = jst::contrib::copyable_box<base_tree_t>;
        using base_node_type = libjst::tree_node_t<base_tree_t>;

        class node_impl;
        class sink_impl;

        wrappee_t _wrappee{};
        std::size_t _max_subtree_depth{};

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr k_depth_tree_impl() = default; //!< Default.

        template <typename wrapped_tree_t>
        explicit constexpr k_depth_tree_impl(wrapped_tree_t && wrappee, std::size_t const max_subtree_depth) noexcept :
            _wrappee{(wrapped_tree_t &&)wrappee},
            _max_subtree_depth{max_subtree_depth}
        {}
        //!\}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(*_wrappee), _max_subtree_depth};
        }

        constexpr sink_impl sink() const noexcept {
            return sink_impl{libjst::sink(*_wrappee)};
        }
   };

    template <typename base_tree_t>
    class k_depth_tree_impl<base_tree_t>::node_impl : public base_node_type {
    private:

        friend k_depth_tree_impl;

        std::size_t _max_subtree_depth{};
        std::size_t _subtree_depth{};

        explicit constexpr node_impl(base_node_type base_node,
                                     std::size_t max_subtree_depth,
                                     std::size_t subtree_depth) noexcept :
            base_node_type{std::move(base_node)},
            _max_subtree_depth{max_subtree_depth},
            _subtree_depth{subtree_depth}
        {}

        explicit constexpr node_impl(base_node_type base_node, std::size_t max_subtree_depth) noexcept :
            node_impl{std::move(base_node), max_subtree_depth, 0u}
        {}

    public:

        explicit constexpr node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const noexcept {
            if (max_depth_reached<true>())
                return std::nullopt;

            return visit<true>(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const noexcept {
            if (max_depth_reached<false>())
                return std::nullopt;

            return visit<false>(base_node_type::next_ref());
        }

    private:

        template <bool is_alt, typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t maybe_child) const {
            if (maybe_child) {
                std::size_t new_depth = _subtree_depth + (on_alternate_path<is_alt>());
                return node_impl{std::move(*maybe_child), _max_subtree_depth, new_depth};
            } else {
                return std::nullopt;
            }
        }

        template <bool is_alt>
        constexpr bool max_depth_reached() const noexcept {
            return (on_alternate_path<is_alt>()) ? _subtree_depth == _max_subtree_depth : false;
        }

        template <bool is_alt>
        constexpr bool on_alternate_path() const noexcept {
            return is_alt || base_node_type::on_alternate_path();
        }

        constexpr friend bool operator==(node_impl const & lhs, sink_impl const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };

    template <typename base_tree_t>
    class k_depth_tree_impl<base_tree_t>::sink_impl {
    private:
        friend k_depth_tree_impl;

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

    namespace _tree_adaptor {
        inline constexpr struct _k_depth
        {
            template <typename base_tree_t, std::unsigned_integral depth_t>
            constexpr auto operator()(base_tree_t && tree, depth_t const depth) const
                noexcept(std::is_nothrow_constructible_v<
                            k_depth_tree_impl<std::remove_reference_t<base_tree_t>>>)
                -> k_depth_tree_impl<std::remove_reference_t<base_tree_t>>
            {
                using adapted_tree_t = k_depth_tree_impl<std::remove_reference_t<base_tree_t>>;
                return adapted_tree_t{(base_tree_t &&)tree, depth};
            }

            template <std::unsigned_integral depth_t>
            constexpr auto operator()(depth_t const depth) const
                noexcept(std::is_nothrow_invocable_v<libjst::tag_t<jst::contrib::make_closure>, depth_t>)
                -> jst::contrib::closure_result_t<_k_depth, depth_t>
            {
                return jst::contrib::make_closure(_k_depth{}, depth);
            }
        } k_depth{};
    } // namespace _tree_adaptor

    using _tree_adaptor::k_depth;
}  // namespace libjst
