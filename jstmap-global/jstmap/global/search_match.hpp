// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides search match.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <optional>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>

#include <jstmap/global/match_position.hpp>

namespace jstmap
{
    struct alignment_result {

        int32_t score{};
        std::vector<seqan3::cigar> cigar_sequence{};

        template <typename original_result_t>
            requires (!std::same_as<std::remove_cvref_t<original_result_t>, alignment_result>)
        explicit alignment_result(original_result_t && res) noexcept :
            score{res.score()},
            cigar_sequence{seqan3::cigar_from_alignment(std::move(res.alignment()))}
        {}
    };

    class search_match
    {
    public:
        search_match() = default;
        explicit search_match(match_position position, alignment_result alignment) noexcept :
            _position{std::move(position)},
            _alignment{std::move(alignment)}
        {}

        void set_position(match_position position) noexcept {
            _position = std::move(position);
        }

        match_position const & position() const noexcept {
            return _position;
        }

        void set_alignment(alignment_result alignment) noexcept {
            _alignment = std::move(alignment);
        }

        bool has_alignment() const noexcept {
            return _alignment.has_value();
        }

        std::vector<seqan3::cigar> get_cigar() const noexcept {
            return _alignment->cigar_sequence;
        }
    private:

        match_position _position{};
        std::optional<alignment_result> _alignment{std::nullopt};
    };
}  // namespace jstmap
