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

#include <jstmap/search/type_alias.hpp>

namespace jstmap
{

    class match {
    private:

        query_index_type _query_id{};

    public:

        match() = default;
        constexpr explicit match(query_index_type const query_id) : _query_id{query_id}
        {}

        constexpr query_index_type query_id() const noexcept {
            return _query_id;
        }
    };
}  // namespace jstmap
