// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides adaption of the field code for the sequence format.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libio/file/field_code.hpp>
#include <libio/format/fasta/fasta_field_code.hpp>
#include <libio/format/fastq/fastq_field_code.hpp>

namespace libio
{
    // define tags for the different fastq fields.
    enum class sequence_field
    {
        // --  0-- forbidden
        id = 1,
        seq,
        qual
    };

    class sequence_field_category
    {
    public:

        constexpr sequence_field_category() = default;
    private:

        template <auto field_tag>
            requires (!std::same_as<decltype(field_tag), sequence_field>)
        constexpr friend bool tag_invoke(tag_t<libio::equivalent>,
                                         sequence_field_category const &,
                                         field_code_type<field_tag> const & fc,
                                         sequence_field const & sequence_ft) noexcept
        {
            switch (sequence_ft)
            {
                case sequence_field::id: return libio::equivalent(fc.category(), fc, fasta_field::id) ||
                                                libio::equivalent(fc.category(), fc, fastq_field::id);
                case sequence_field::seq: return libio::equivalent(fc.category(), fc, fasta_field::seq) ||
                                                 libio::equivalent(fc.category(), fc, fastq_field::seq);
                case sequence_field::qual: return libio::equivalent(fc.category(), fc, fastq_field::qual);
                default: return false; // no known equivalency!
            }
        }
    };

    template <>
    struct field_code_category<sequence_field> : std::type_identity<sequence_field_category>
    {
    };

    static_assert(libio::equivalent(default_field_code_category{}, field_code<fasta_field::id>, fasta_field::id));

} // namespace libio
