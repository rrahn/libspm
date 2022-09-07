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

#include <seqan3/core/debug_stream.hpp>

#include <libjst/matcher/horspool_matcher.hpp>
#include <libjst/matcher/shiftor_matcher.hpp>

#include <libjst/search/polymorphic_sequence_searcher.hpp>
#include <libjst/referentially_compressed_sequence_store/debug_rcs_store.hpp>

// #include <libjst/search/concept.hpp>
// #include <libjst/search/search_base.hpp>
// #include <libjst/traversal/searcher_factory.hpp>
// #include <libjst/search/shift_or_search.hpp>
// #include <libjst/search/horspool_search.hpp>
// #include <libjst/search/myers_search.hpp>
// #include <libjst/search/pigeonhole_filter.hpp>
// #include <libjst/search/state_manager_stack.hpp>

#include <jstmap/search/search_queries.hpp>

namespace jstmap
{

std::vector<search_match2> search_queries_horspool(rcs_store_t const & jst, bin_t const & queries, float const error_rate)
{
    assert(!seqan::empty(queries));


    std::vector<search_match2> matches{};
    size_t query_id{};
    std::ranges::for_each(queries, [&] (auto const & query)
    {
        libjst::horspool_matcher matcher{query};

        libjst::polymorphic_sequence_searcher searcher{jst};
        searcher(matcher, [&] (auto && finder, [[maybe_unused]] auto && jst_context) {
            // auto && [node, finder] = hit; // how can we now check the availability of the position?
            matches.emplace_back(seqan::beginPosition(finder), seqan::endPosition(finder), query_id, 0);
        });
        ++query_id;
    });
    return matches;
}

std::vector<search_match2> search_queries_shiftor(rcs_store_t const & jst, bin_t const & queries, float const error_rate)
{
    assert(!seqan::empty(queries));

    // libjst::debug_rcs_store _jst{jst, {0, jst.source().size()}, {0, jst.variants().size()}};

    std::vector<search_match2> matches{};
    size_t query_id{};
    // std::ranges::for_each(queries, [&] (auto const & query)
    for (auto const & query : queries)
    {
        // auto const & query = queries[1];
        libjst::shiftor_matcher matcher{query};

        auto callback = [&] (auto && finder, [[maybe_unused]] auto && jst_context) {
            // auto && [node, finder] = hit; // how can we now check the availability of the position?
            matches.emplace_back(seqan::beginPosition(finder), seqan::endPosition(finder), query_id, 0);

        };

        libjst::polymorphic_sequence_searcher searcher{jst};
        searcher(matcher, callback);

        ++query_id;
    }
    //);
    return matches;
}

}  // namespace jstmap
