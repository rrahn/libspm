// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/copyable_box.hpp>

#include <libjst/sequence_tree/concept.hpp>

namespace libjst
{
    template <typename base_tree_t,
              template <typename, typename, typename...> typename node_extension_t,
              typename ...extension_args_t>
    class extendable_tree {
    private:

        using base_node_type = libjst::tree_node_t<base_tree_t>;
        using sink_type = libjst::tree_sink_t<base_tree_t>;

        class node_impl;

        using tree_box_t = jst::contrib::copyable_box<base_tree_t>;

        base_tree_t _wrappee{};

    public:

        template <typename wrappee_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, extendable_tree> &&
                      std::constructible_from<base_tree_t, wrappee_t>)
        constexpr explicit extendable_tree(wrappee_t && wrappee) noexcept : _wrappee{(wrappee_t &&)wrappee}
        {}

        constexpr node_impl root() const noexcept {
            return node_impl{libjst::root(_wrappee)};
        }

        constexpr sink_type sink() const noexcept {
            return libjst::sink(_wrappee);
        }
    };

    template <typename base_tree_t,
              template <typename, typename, typename...> typename node_extension_t,
              typename ...extension_args_t>
    class extendable_tree<base_tree_t, node_extension_t, extension_args_t...>::node_impl :
            public base_node_type,
            public node_extension_t<node_impl, base_node_type, extension_args_t...> {
    private:
        using extension_t = node_extension_t<node_impl, base_node_type, extension_args_t...>;

        friend extension_t;
        friend extendable_tree;

        explicit constexpr node_impl(base_node_type && base_node) noexcept :
            base_node_type{std::move(base_node)}
        {
            extension_t::initialise();
        }

       explicit constexpr node_impl(base_node_type && base_node, extension_t extension) noexcept :
            base_node_type{std::move(base_node)},
            extension_t{std::move(extension)}
        {}

    public:

        node_impl() = default;

        constexpr std::optional<node_impl> next_alt() const {
            return visit(base_node_type::next_alt());
        }

        constexpr std::optional<node_impl> next_ref() const {
            return visit(base_node_type::next_ref());
        }

        constexpr auto operator*() const noexcept {
            return extension_t::label(*static_cast<base_node_type const &>(*this));
        }

    private:

        template <typename maybe_child_t>
        constexpr std::optional<node_impl> visit(maybe_child_t && maybe_child) const {
            if (maybe_child) {
                extension_t child_extension = extension_t::notify(*maybe_child);
                return node_impl{std::move(*maybe_child), std::move(child_extension)};
            } else {
                return std::nullopt;
            }
        }

        friend bool operator==(node_impl const & lhs, sink_type const & rhs) noexcept {
            return static_cast<base_node_type const &>(lhs) == rhs;
        }
    };
}  // namespace libjst
