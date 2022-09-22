// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the stream token.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libio/file/tokenization.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    class noop_token
    {
    public:

        noop_token() = default;

    private:

        template <typename value_t>
        constexpr friend auto tag_invoke(tag_t<libio::detokenize_to>, noop_token const &, value_t const &) noexcept
        { // noop
        }
    };
} // namespace libio
