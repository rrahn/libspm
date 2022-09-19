// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the sequence_record.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <type_traits>

#include <libio/format/sequence/sequence_field_code.hpp>
#include <libio/file/tokenization.hpp>

namespace libio
{
    class sequence_record
    {
    private:
        std::string _id{};
        std::string _seq{};
        std::string _qual{};

    public:
        constexpr sequence_record() = default;

        std::string const & id() const noexcept
        {
            return _id;
        }

        std::string const & seq() const noexcept
        {
            return _seq;
        }

        std::string const & qual() const noexcept
        {
            return _qual;
        }

    private:

        using category_t = field_code_category_t<sequence_field>;

        template <auto field_tag, typename chunked_buffer_t>
            requires (libio::equivalent(category_t{}, field_code<field_tag>, sequence_field::id))
        friend auto tag_invoke(tag_t<libio::set_field>,
                               sequence_record &me,
                               field_code_type<field_tag>,
                               chunked_buffer_t &&buffer)
        {
            libio::read_token(me._id, buffer);
        }

        template <auto field_tag, typename chunked_buffer_t>
            requires (libio::equivalent(category_t{}, field_code<field_tag>, sequence_field::seq))
        friend auto tag_invoke(tag_t<libio::set_field>,
                               sequence_record &me,
                               field_code_type<field_tag>,
                               chunked_buffer_t &&buffer)
        {
            libio::read_token(me._seq, buffer);
        }

        template <auto field_tag, typename chunked_buffer_t>
            requires (libio::equivalent(category_t{}, field_code<field_tag>, sequence_field::qual))
        friend auto tag_invoke(tag_t<libio::set_field>, sequence_record &me,
                               field_code_type<field_tag>,
                               chunked_buffer_t &&buffer)
        {
            libio::read_token(me._qual, buffer);
        }
    };
} // namespace libio
