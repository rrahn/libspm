// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the fasta format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>
#include <vector>
#include <string>

#include <seqan/seq_io.h>

#include <libio/format/fasta/fasta_record.hpp>
#include <libio/format/format_concept.hpp>
#include <libio/format/format_extension.hpp>
#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    class fasta_format : public format_extension
    {
    public:

        fasta_format() : format_extension{".fa", ".fasta", ".fn"}
        {
        }

    private:

        template <typename stream_t>
        friend fasta_record tag_invoke(tag_t<libio::get_record>, fasta_format const &/*format*/, stream_t & stream)
        {
            std::string _seq{};
            std::string _id{};
            seqan::readRecord(_id, _seq, stream, seqan::Fasta{}); // read with tokenisation!
            return libio::fasta_record{std::move(_id), std::move(_seq)};
        }

        friend seqan::Fastq tag_invoke(tag_t<libio::format_tag>, fasta_format const & /*format*/) noexcept
        {
            return seqan::Fastq{};
        }
    };
} // namespace libio
