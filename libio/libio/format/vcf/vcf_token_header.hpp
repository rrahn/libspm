// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementations of a combined token.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libio/file/tokenization.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    template <typename meta_token_t, typename sample_token_t>
    class vcf_token_header
    {
    private:
        meta_token_t _meta_token{};
        sample_token_t _sample_token{};

    public:
        vcf_token_header() = delete;
        vcf_token_header(meta_token_t && meta_token, sample_token_t && sample_token) noexcept :
            _meta_token{std::move(meta_token)},
            _sample_token{std::move(sample_token)}
        {}

    private:

        template <typename header_t>
        friend void tag_invoke(tag_t<libio::detokenize_to>, vcf_token_header &me, header_t &header)
        {
            { // make sure first token is consumed before second token.
                meta_token_t token = std::move(me._meta_token);
                libio::detokenize_to(token, header);
            }

            {
                sample_token_t token = std::move(me._sample_token);
                libio::detokenize_to(token, header);
            }
        }
    };
}  // namespace libio
