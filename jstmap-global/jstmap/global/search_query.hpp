// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides query type.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <jstmap/global/jstmap_types.hpp>

namespace jstmap
{

    class search_query {
    public:

        using key_type = std::size_t;
        using value_type = sequence_record_t;

        search_query() = default;
        search_query(key_type const key, value_type value)
            noexcept(std::is_nothrow_move_constructible_v<value_type>) : _key{key}, _value{std::move(value)}
        {}

        key_type key() const noexcept {
            return _key;
        }
        value_type value() const noexcept {
            return _value;
        }
    private:

        key_type _key{};
        value_type _value{};
    };
}  // namespace jstmap
