// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implements find_if functionality.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <ranges>

namespace libio
{

    template <std::ranges::random_access_range range_t, typename predicate_t>
        requires std::predicate<predicate_t, std::ranges::range_reference_t<range_t>>
    inline std::ranges::iterator_t<range_t> find_if(range_t &&range, predicate_t &&pred) noexcept
    {
        for (auto it = std::ranges::begin(range); it != std::ranges::end(range); ++it)
        {
            if (pred(*it))
                return it;
        }

        return std::ranges::end(range);
    }

    // template <std::ranges::random_access_range range_t, typename predicate_t>
    //     requires std::predicate<predicate_t, std::ranges::range_reference_t<range_t>>
    // inline std::ranges::iterator_t<range_t> find_if_not(range_t && range, predicate_t && pred) noexcept
    // {
    //     int steps = 0;
    //     for (; steps < std::ranges::ssize(range); ++steps)
    //     {
    //         if (!pred(range[steps]))
    //             break;
    //     }

    //     return std::ranges::next(std::ranges::begin(range), steps);
    // }
} // namespace libio
