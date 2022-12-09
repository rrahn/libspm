// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides node traits class.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <typename node_t>
    struct rcs_node_traits
    {};

    template <typename node_t>
        requires requires { typename node_t::rcs_store_type; }
    struct rcs_node_traits<node_t> {
        using rcs_store_type = typename node_t::rcs_store_type;
        using sequence_type = typename rcs_store_type::source_type;
        using variant_map_type = typename rcs_store_type::variant_map_type;
        using variant_iterator = std::ranges::iterator_t<variant_map_type const>;
        using variant_type = std::ranges::range_value_t<variant_map_type const>;
        using breakpoint_type = libjst::variant_breakpoint_t<variant_type>;
    };
}  // namespace libjst
