// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for the sequence tree classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Operation CPOs for sequence trees
    // ----------------------------------------------------------------------------

    // root
    namespace _root {
        inline constexpr struct _cpo  {
            template <typename tree_t>
                requires std::tag_invocable<_cpo, tree_t>
            constexpr auto operator()(tree_t && tree) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, tree_t>)
                -> std::tag_invoke_result_t<_cpo, tree_t>
            {
                return std::tag_invoke(_cpo{}, (tree_t &&) tree);
            }

        private:

            template <typename tree_t>
                requires requires (tree_t tree) { { std::forward<tree_t &&>(tree).root() }; }
            constexpr friend auto tag_invoke(_cpo, tree_t && tree)
                noexcept(noexcept(std::declval<tree_t>().root()))
                -> decltype(std::declval<tree_t>().root())
            {
                return ((tree_t &&)tree).root();
            }
        } root;
    } // namespace _root
    using _root::root;

    template <typename tree_t>
    using tree_node_t = std::invoke_result_t<_root::_cpo, tree_t>;

    // sink
    namespace _sink {
        inline constexpr struct _cpo  {
            template <typename tree_t>
                requires std::tag_invocable<_cpo, tree_t>
            constexpr auto operator()(tree_t && tree) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, tree_t>)
                -> std::tag_invoke_result_t<_cpo, tree_t>
            {
                return std::tag_invoke(_cpo{}, (tree_t &&)tree);
            }
        private:

            template <typename tree_t>
                requires requires (tree_t tree) { { std::forward<tree_t &&>(tree).sink() }; } // has member function sink!
            constexpr friend auto tag_invoke(_cpo, tree_t && tree)
                noexcept(noexcept(std::declval<tree_t>().sink()))
                -> decltype(std::declval<tree_t>().sink())
            {
                return ((tree_t &&)tree).sink();
            }
        } sink;
    } // namespace _sink
    using _sink::sink;

    template <typename tree_t>
    using tree_sink_t = std::invoke_result_t<_sink::_cpo, tree_t>;

    // ----------------------------------------------------------------------------
    // Label type of trees and tree nodes.
    // ----------------------------------------------------------------------------
    template <typename node_t>
    using node_label_t = decltype(*std::declval<node_t>());

    template <typename tree_t>
    using tree_label_t = node_label_t<tree_node_t<tree_t>>;

    // ----------------------------------------------------------------------------
    // Operation CPOs for sequence tree nodes
    // ----------------------------------------------------------------------------

    namespace _next_alt {
        inline constexpr struct _cpo  {
            template <typename node_t>
                requires std::tag_invocable<_cpo, node_t>
            constexpr auto operator()(node_t && node) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, node_t>)
                -> std::tag_invoke_result_t<_cpo, node_t>
            {
                return std::tag_invoke(_cpo{}, (node_t &&)node);
            }

        private:

            // has member function next_alt!
            template <typename node_t>
                requires requires (node_t && node) { { std::forward<node_t &&>(node).next_alt() }; }
            constexpr friend auto tag_invoke(_cpo, node_t && node)
                noexcept(noexcept(std::declval<node_t &&>().next_alt()))
                -> decltype(std::declval<node_t &&>().next_alt())
            {
                return std::forward<node_t &&>(node).next_alt();
            }
        } next_alt;
    } // namespace _next_alt
    using _next_alt::next_alt;

    namespace _next_ref {
        inline constexpr struct _cpo  {
            template <typename node_t>
                requires std::tag_invocable<_cpo, node_t>
            constexpr auto operator()(node_t && node) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, node_t>)
                -> std::tag_invoke_result_t<_cpo, node_t>
            {
                return std::tag_invoke(_cpo{}, (node_t &&)node);
            }

        private:

            // has member function next_ref!
            template <typename node_t>
                requires requires (node_t && node) { { std::forward<node_t &&>(node).next_ref() }; }
            constexpr friend auto tag_invoke(_cpo, node_t && node)
                noexcept(noexcept(std::declval<node_t &&>().next_ref()))
                -> decltype(std::declval<node_t &&>().next_ref())
            {
                return std::forward<node_t &&>(node).next_ref();
            }
        } next_ref;
    } // namespace _next_ref
    using _next_ref::next_ref;

    namespace _prev_alt {
        inline constexpr struct _cpo  {
            template <typename node_t>
                requires std::tag_invocable<_cpo, node_t>
            constexpr auto operator()(node_t && node) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, node_t>)
                -> std::tag_invoke_result_t<_cpo, node_t>
            {
                return std::tag_invoke(_cpo{}, (node_t &&)node);
            }

        private:

            // has member function prev_alt!
            template <typename node_t>
                requires requires (node_t && node) { { std::forward<node_t &&>(node).prev_alt() }; }
            constexpr friend auto tag_invoke(_cpo, node_t && node)
                noexcept(noexcept(std::declval<node_t &&>().prev_alt()))
                -> decltype(std::declval<node_t &&>().prev_alt())
            {
                return std::forward<node_t &&>(node).prev_alt();
            }
        } prev_alt;
    } // namespace _prev_alt
    using _prev_alt::prev_alt;

    namespace _prev_ref {
        inline constexpr struct _cpo  {
            template <typename node_t>
                requires std::tag_invocable<_cpo, node_t>
            constexpr auto operator()(node_t && node) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, node_t>)
                -> std::tag_invoke_result_t<_cpo, node_t>
            {
                return std::tag_invoke(_cpo{}, (node_t &&)node);
            }

        private:

            // has member function prev_ref!
            template <typename node_t>
                requires requires (node_t && node) { { std::forward<node_t &&>(node).prev_ref() }; }
            constexpr friend auto tag_invoke(_cpo, node_t && node)
                noexcept(noexcept(std::declval<node_t &&>().prev_ref()))
                -> decltype(std::declval<node_t &&>().prev_ref())
            {
                return std::forward<node_t &&>(node).prev_ref();
            }
        } prev_ref;
    } // namespace _prev_ref
    using _prev_ref::prev_ref;

    // ----------------------------------------------------------------------------
    // Tree concepts
    // ----------------------------------------------------------------------------

    // namespace detail::exposition_only {

    //     template <typename bool_t>
    //     concept boolean_testable_impl = std::convertible_to<bool_t, bool>;

    //     template <typename bool_t>
    //     concept boolean = boolean_testable_impl<bool_t> &&
    //                       requires (bool_t && b) { { !std::forward<bool_t>(b) } -> boolean_testable_impl; };

    //     template<typename t1, typename t2>
    //     concept weakly_equality_comparable_with =
    //         requires (std::remove_reference_t<t1> const & lhs, std::remove_reference_t<t2> const & rhs) {
    //         { lhs == rhs } -> boolean;
    //         { lhs != rhs } -> boolean;
    //         { rhs == lhs } -> boolean;
    //         { rhs != lhs } -> boolean;
    //     };

    //     template <typename nil_t, typename node_t>
    //     concept nullable_node_for = requires (nil_t node) {
    //                                     { !node } -> boolean;
    //                                     { *node } -> std::convertible_to<node_t>;
    //                                 };
    // } // namespace detail::exposition_only

    // template <typename node_t>
    // concept directed_node = std::semiregular<node_t> &&
    //                         requires (node_t n) {
    //                             { libjst::next_alt(std::forward<node_t>(n)) } -> detail::exposition_only::nullable_node_for<node_t>;
    //                             { libjst::next_ref(std::forward<node_t>(n)) } -> detail::exposition_only::nullable_node_for<node_t>;
    //                             // { libjst::is_leaf(std::as_const(n)) } -> detail::exposition_only::boolean;
    //                         };


    // template <typename node_t>
    // concept undirected_node = directed_node<node_t> &&
    //                           requires (node_t n) {
    //                             { libjst::prev_alt(n) } -> detail::exposition_only::nullable_node_for<node_t>;
    //                             { libjst::prev_ref(n) } -> detail::exposition_only::nullable_node_for<node_t>;
    //                           };

    // template <typename sink_t, typename node_t>
    // concept tree_sink_for = std::semiregular<sink_t> &&
    //                         directed_tree_node<node_t> &&
    //                         detail::exposition_only::weakly_equality_comparable_with<sink_t, node_t>;


    // template <typename tree_t>
    // concept rooted_tree = requires (tree_t const & t) {
    //     std::same_as
    //     { libjst::root(t) } -> directed_tree_node;
    //     { libjst::sink(t) } -> directed_tree_node;
    // };

    // template <typename tree_t>
    // concept sequence_tree = rooted_tree<tree_t> &&
    //                         labeled_tree<tree_t> &&
    //                         requires (tree_t const & t) {
    //                             typename libjst::node_label_t<tree_t>;
    //                         }

    // template <typename tree_t>
    // concept directed_sequence_tree = sequence_tree<tree_t> && directed_node<libjst::tree_node_t<tree_t>>;

    // template <typename tree_t>
    // concept undirected_sequence_tree = directed_sequence_tree<tree_t> && undirected_node<libjst::tree_node_t<tree_t>>;

}  // namespace libjst
