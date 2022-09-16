// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the fasta_record.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <type_traits>

#include <libio/format/fastq/fastq_field_code.hpp>
#include <libio/file/tokenization.hpp>

namespace libio
{
    // Implementation of some custom record type.
    class fastq_record
    {
    private:
        std::string _value{};
        size_t _seq_offset{};
        size_t _qual_offset{};

    public:

        std::string_view seq() const noexcept
        {
            return std::string_view{_value.data() + _seq_offset, _qual_offset - _seq_offset};
        }

        std::string_view id() const noexcept
        {
            return std::string_view{_value.data(), _seq_offset};
        }

        std::string_view qual() const noexcept
        {
            return std::string_view{_value.data() + _qual_offset, _value.size() - _qual_offset};
        }

    private:
        template <typename chunked_buffer_t>
        friend auto tag_invoke(tag_t<libio::set_field>, fastq_record & me,
                                                        field_code_type<fastq_field::id>,
                                                        chunked_buffer_t && buffer) {
            libio::read_token(me._value, buffer); //fill record data
            me._seq_offset = me._value.size();
        }

        template <typename chunked_buffer_t>
        friend auto tag_invoke(tag_t<libio::set_field>, fastq_record & me,
                                                        field_code_type<fastq_field::seq>,
                                                        chunked_buffer_t && buffer) {
            libio::read_token(me._value, buffer);
            me._qual_offset = me._value.size(); // dump how much was written into the buffer.
        }

        template <typename chunked_buffer_t>
        friend auto tag_invoke(tag_t<libio::set_field>, fastq_record & me,
                                                        field_code_type<fastq_field::qual>,
                                                        chunked_buffer_t && buffer) {
            libio::read_token(me._value, buffer);
        }
    };
} // namespace libio
