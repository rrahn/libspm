// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <iosfwd>
#include <memory>

#include <libio/formats/sequence/sequence_record_raw.hpp>

namespace libio
{
    class sequence_format
    {
    private:
        struct format_base;

        std::unique_ptr<format_base> _format{};

    public:
        // we are called with a specific list of formats
        template <typename format_t>
        // requires sequence_format<format_t>?
        sequence_format(format_t &&format) noexcept(std::is_nothrow_constructible_v<format_impl<format_t>, format_t>)
            : _format{std::make_unique<format_impl<format_t>>(std::forward<format_t>(format))}
        {
        }
        // fasta_format, fastq_format, genbank_format, sam_sequence_format, bam_sequence_format ...
        // maybe not sure which one it is?
        // now some one calls read_record
        // which format is set? and when is it set?
        // Doesn't implement read/write header because it doesn't have one.

        sequence_record_raw read_record(std::istream &istream)
        {
            return _format->read_record(istream);
        }

        template <typename record_t>
        // requires sequence_record<record_t>
        void write_record(std::ostream &ostream, record_t &&record)
        {
            // converting via CPO to a record: which night even be type erased in here.
            // return _format->write_record(ostream);
        }

    private:
        // type erased format
        struct format_base
        {

            virtual sequence_record_raw read_record(std::istream &) { return sequence_record_raw{}; };
            // TODO: write record
        };

        template <typename format_t>
        struct format_impl final : public format_base
        {
            format_t _format{};

            format_impl() = delete;
            format_impl(format_t format) : _format{std::forward<format_t>(format)}
            {
            }

            sequence_record_raw read_record(std::istream &istream) override
            {
                return _format.read_record(istream);
            }
        };
    };
} // namespace libio
