// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides alternate sequence kind enumeration.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <string_view>
#include <type_traits>
#include <iosfwd>

namespace libjst
{
    enum class alternate_sequence_kind {
        insertion = 0,
        replacement = 1,
        deletion = 2
    };

    template <typename char_t, typename char_traits_t = std::char_traits<char_t>>
    inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                                alternate_sequence_kind const & alt_kind) {

        auto to_message = [] (alternate_sequence_kind const & kind) {
            using namespace std::literals;
            switch (kind) {
                case alternate_sequence_kind::replacement: return "replacement"sv;
                case alternate_sequence_kind::deletion: return "deletion"sv;
                case alternate_sequence_kind::insertion: return "insertion"sv;
                default: return "unknown"sv;
            };
        };
        return stream << to_message(alt_kind);
    }
}  // namespace libjst
