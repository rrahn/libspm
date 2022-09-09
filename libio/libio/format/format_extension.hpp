// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the extension class.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>
#include <string>

namespace libio
{
    class format_extension
    {
    private:
        std::vector<std::string> _value;

    public:
        format_extension() = delete;
        format_extension(auto &&...extensions) : _value{extensions...}
        {
        }

        std::vector<std::string> const &valid_extensions() const noexcept
        {
            return _value;
        }
    };
} // namespace libio
