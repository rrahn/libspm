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

std::vector<libjst::context_position> search_queries(partitioned_jst_t const & partitioned_jst, std::vector<raw_sequence_t> const & queries)
{
    assert(!queries.empty());

    std::vector<libjst::context_position> results{};

    std::ranges::for_each(queries, [&] (raw_sequence_t const & query)
    {
        // prepare searcher
        using state_t = typename decltype(libjst::shift_or_pattern_searcher{query})::state_type;
        libjst::shift_or_pattern_searcher searcher{query, libjst::search_state_manager_stack<state_t>{}};

        for (uint32_t index = 0; index < partitioned_jst.bin_count(); ++index)
        {
            auto jst_range_agent = partitioned_jst.range_agent(libjst::context_size{static_cast<uint32_t>(query.size())},
                                                               libjst::bin_index{index},
                                                               searcher.state_manager()); // already pushing a branch.
            searcher(jst_range_agent, [&] (auto & it)
            {
                process_hits(results, jst.sequence_positions_at(it.coordinate()));
            });
        }
    });

    return results;
}

}  // namespace jstmap
