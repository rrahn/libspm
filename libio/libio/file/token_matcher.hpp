// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of a token matcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <string_view>

namespace libio
{

    class token_matcher
    {
    private:
        std::string_view _pattern{};


    public:
        token_matcher() = default;
        token_matcher(std::string_view pattern) : _pattern{std::move(pattern)}
        {
        }

        // constexpr std:: operator()(std::string_view text) const noexcept // what do we return?
        // {
        //     // last state?
        //     // number of hits?

        // }
    };

} // namespace libio
