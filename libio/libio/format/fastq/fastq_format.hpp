// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the fastq format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>

#include <libio/format/fastq/fastq_token.hpp>
#include <libio/format/format_concept.hpp>
#include <libio/format/format_extension.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    class fastq_format : public format_extension
    {
    public:

        fastq_format() : format_extension{".fq", ".fastq"}
        {
        }

    private:

        template <typename char_t, typename char_traits_t>
        constexpr friend auto tag_invoke(tag_t<libio::format_token>,
                                         fastq_format const &,
                                         std::basic_istream<char_t, char_traits_t> & istream)
        {
            return fastq_token{istream};
        }
    };
} // namespace libio
