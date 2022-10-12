// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides build function to create the journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <jstmap/search/type_alias.hpp>

#include <libjst/traversal/jst_node.hpp>

namespace jstmap
{

using jst_bin_t = typename partitioned_jst_t::traverser_model_t;

struct search_match
{
    using journal_decorator_type = typename jst_t::journal_decorator_type;
    using journal_decorator_iterator_type = std::ranges::iterator_t<journal_decorator_type>;
    using subrange_type = std::ranges::subrange<journal_decorator_iterator_type, journal_decorator_iterator_type>;
    // query_id; -> ordered vector by query_id in order to filter duplicates
    // jst_span; -> sequence region identified within the jst
    // error_count; -> number of errors identified
    // jst_coordinate; -> base coordinate where the hit was found (is this enough for overlap duplicates?)
    journal_decorator_type jst_sequence{};
    size_t begin_position{};
    size_t end_position{};
    libjst::journal_sequence_tree_coordinate hit_coordinate{};
    size_t query_id{};
    size_t error_count{};

    subrange_type sequence() const noexcept
    {
        return subrange_type{jst_sequence.begin() + begin_position, jst_sequence.begin() + end_position};
    }

    auto operator<=>(search_match const & rhs) const noexcept
    {
        return std::tie(query_id, error_count) <=> std::tie(rhs.query_id, rhs.error_count);
    }
};

struct search_match2
{
    using node_t = libjst::jst_node<fwd_jst_t>;

    // jst_span; -> sequence region identified within the jst
    // error_count; -> number of errors identified
    // jst_coordinate; -> base coordinate where the hit was found (is this enough for overlap duplicates?)
    node_t node{}; // do we copy the node?
    size_t first_position{};
    size_t last_position{};
    size_t query_id{};
    size_t error_count{};

    auto sequence() const noexcept
    {
        auto const & seq = node.sequence();
        return std::ranges::subrange{seq.begin() + first_position, seq.begin() + last_position};
    }

    auto operator<=>(search_match2 const & rhs) const noexcept
    {
        return std::tie(query_id, error_count) <=> std::tie(rhs.query_id, rhs.error_count);
    }
};

std::vector<search_match> search_queries_(jst_bin_t const &, bin_t const &, float const);
std::vector<libjst::context_position> search_queries(partitioned_jst_t const &, std::vector<raw_sequence_t> const &);

std::vector<search_match2> search_queries_horsppol(fwd_jst_t const &, bin_t const &, float const);
std::vector<search_match2> search_queries_shiftor(fwd_jst_t const &, bin_t const &, float const);

}  // namespace jstmap
