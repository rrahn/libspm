// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <memory>
#include <string_view>

namespace libio
{
    class sequence_record_raw
    {
    public:
        using seq_t = std::string_view;
        using id_t = std::string_view;
        using qual_t = std::string_view;

        sequence_record_raw() = default;
        template <typename record_t>
        // requires sequence_record
        sequence_record_raw(record_t &&sequence_record) noexcept(std::is_nothrow_constructible_v<record_t, record_t &&>)
            : _record{std::make_unique<record_impl<record_t>>(std::forward<record_t>(sequence_record))}
        {
        }

        seq_t seq() const
        {
            return _record->seq();
        }

        id_t id() const
        {
            return _record->id();
        }

        qual_t qual() const
        {
            return _record->qual();
        }

    private:
        // type erased record
        struct record_base
        {
            virtual seq_t seq() const { return seq_t{}; };
            virtual id_t id() const { return id_t{}; };
            virtual qual_t qual() const { return qual_t{}; };
        };

        template <typename record_t>
        struct record_impl final : public record_base
        {
            record_t _record{};

            record_impl() = delete;
            record_impl(record_t record) : _record{std::forward<record_t>(record)}
            {
            }

            seq_t seq() const override
            {
                return _record.seq();
            }

            id_t id() const override
            {
                return _record.id();
            }

            qual_t qual() const override
            {
                return _record.qual();
            }
        };

        std::unique_ptr<record_base> _record{};
    };
} // namespace libio
