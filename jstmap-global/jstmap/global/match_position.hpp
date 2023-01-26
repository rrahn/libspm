// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a reference position.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iosfwd>

#include <libjst/sequence_tree/seek_position.hpp>

namespace jstmap
{
    struct match_position {
    public:
        libjst::seek_position tree_position{};
        std::ptrdiff_t label_offset{};

    private:

        constexpr friend bool operator==(match_position const &, match_position const &) noexcept = default;
        constexpr friend auto operator<=>(match_position const &, match_position const &) noexcept = default;
    };

    template <typename char_t, typename char_traits_t, typename match_position_t>
        requires std::same_as<std::remove_cvref_t<match_position_t>, match_position>
    inline std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                                  match_position_t && position)
    {
        stream << "<tree_position = " << position.tree_position << " label_offset = " << position.label_offset << ">";
        return stream;
    }
}  // namespace jstmap
