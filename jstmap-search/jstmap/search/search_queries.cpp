// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <iterator>
#include <tuple>
#include <ranges>

#include <libjst/search/shift_or_search.hpp>
#include <libjst/search/state_manager_stack.hpp>
#include <libjst/journaled_sequence_tree.hpp>

#include <jstmap/search/write_results.hpp>
#include <jstmap/search/search_queries.hpp>

namespace jstmap
{

void process_hits(std::vector<libjst::context_position> & results,
                  std::vector<libjst::context_position> const & cursor_positions)
{
    std::ranges::for_each(cursor_positions, [&] (libjst::context_position const & mapping_position)
    {
        results.push_back(mapping_position);
    });
}

std::vector<libjst::context_position> search_queries(jst_t const & jst, std::vector<raw_sequence_t> const & queries)
{
    assert(!queries.empty());

    std::vector<libjst::context_position> results{};

    std::ranges::for_each(queries, [&] (raw_sequence_t const & query)
    {
        using state_t = typename decltype(libjst::shift_or_pattern_searcher{query})::state_type;
        using state_manager_t = libjst::search_state_manager_stack<state_t>;
        libjst::shift_or_pattern_searcher searcher{query, state_manager_t{}};

        auto range_agent = jst.range_agent(query.size(), searcher.state_manager());
        searcher(range_agent, [&] (auto & it) { process_hits(results, it.positions()); });
    });

    return results;
}

}  // namespace jstmap
