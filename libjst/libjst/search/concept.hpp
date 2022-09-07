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

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Operation CPOs
    // ----------------------------------------------------------------------------



    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    // template <typename jst_t>
    // concept journaled_sequence_tree = requires (jst_t const & jst)
    // {
    //     typename libjst_t::base_sequence_t<jst_t>;
    //     typename libjst_t::variant_store_t<jst_t>;
    //     // TODO: common type between reference and variant sequence

    //     { libjst::variant_store(jst) } -> libjst::covered_sequence_variant_store;
    //     { libjst::base_sequence(jst) } -> std::ranges::random_access_range;
    // };

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    // template <typename searcher_t, typename text_t>
    // concept searcher_for = requires
    // (searcher_t const & op, text_t const & text)
    // {
    //     { libjst::window_size(op) } -> std::integral;
    //     { op(text, [] (auto &&) { }) };
    // };

    // template <typename searcher_t, typename text_t>
    // concept stateful_searcher_for = searcher_for<searcher_t, text_t> && requires
    // (searcher_t const & c_op, searcher_t & op)
    // {
    //     typename libjst::search_state_t<searcher_t>;
    //     { libjst::state(c_op) } -> std::same_as<libjst::search_state_t<searcher_t> const &>; // TODO: must be searcher for
    //     { libjst::state(op) } -> std::same_as<libjst::search_state_t<searcher_t> &>; // TODO: must be searcher for
    //     { libjst::state(op, libjst::search_state_t<searcher_t>{}) };
    // };
}  // namespace libjst
