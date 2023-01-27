// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides match.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <jstmap/global/search_match.hpp>
#include <jstmap/global/search_query.hpp>

namespace jstmap
{

    class search_matches {
    public:

        using search_query_type = jstmap::search_query;
        using match_type = search_match;
        using match_list_type = std::vector<match_type>;

        search_matches() = default;
        explicit search_matches(search_query_type query) noexcept(std::is_nothrow_move_constructible_v<search_query_type>) :
            _query{std::move(query)}
        {}

        search_query_type const & query() const noexcept {
            return _query;
        }

        match_list_type const & matches() const noexcept {
            return _match_list;
        }

        void record_match(match_type match) noexcept {
            _match_list.push_back(std::move(match));
        }

    private:

        search_query_type _query{};
        match_list_type _match_list{};
    };
}  // namespace jstmap
