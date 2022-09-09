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

namespace libio
{
    class fasta_record
    {
    private:
        std::string _value{};
        size_t _seq_offset{};

    public:
        fasta_record() = default;
        // TODO: We need a special value type here as well.
        fasta_record(std::string value) noexcept(std::is_nothrow_move_constructible_v<std::string>)
            : _value{std::move(value)}
        {
        }

        fasta_record(std::string id, std::string seq) noexcept(std::is_nothrow_move_assignable_v<std::string>)
        {
            _seq_offset = id.size(),
            _value = std::move(id) + std::move(seq);
        }

        std::string_view seq() const noexcept
        {
            return std::string_view{_value.data() + _seq_offset, _value.size() - _seq_offset};
        }

        std::string_view id() const noexcept
        {
            return std::string_view{_value.data(), _seq_offset};
        }
    };
} // namespace libio
