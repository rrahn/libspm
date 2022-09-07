// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides CPOs for the search interface.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename node_label_t>
    concept journaled_sequence_tree_node_label = std::semiregular<std::remove_cvref_t<node_label_t>> &&
    requires (node_label_t & label)
    {
        { std::as_const(label).sequence() };
        { std::as_const(label).coverage() };
        { std::as_const(label).has_value() } -> std::same_as<bool>;
        { !std::as_const(label) } -> std::same_as<bool>;
    };

    template <typename node_t>
    concept journaled_sequence_tree_node = std::semiregular<std::remove_cvref_t<node_t>> && requires (node_t & node)
    {
        typename std::remove_cvref_t<node_t>::label_type;

        { *std::as_const(node) } -> journaled_sequence_tree_node_label;
        { std::as_const(node).is_leaf() } -> std::same_as<bool>;
        { !std::as_const(node) } -> std::same_as<bool>;
        { node.branch() } -> std::same_as<bool>;
    };
}  // namespace libjst
