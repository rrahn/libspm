// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of specific snp sequence variant.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <span>

#include <seqan3/alphabet/concept.hpp>

#include <libjst/sequence_variant/concept.hpp>

namespace libjst
{
    template <seqan3::semialphabet alphabet_t>
    inline const std::array<alphabet_t, seqan3::alphabet_size<alphabet_t>> snp_value{[] ()
    {
        std::array<alphabet_t, seqan3::alphabet_size<alphabet_t>> tmp{};
        using rank_t = seqan3::alphabet_rank_t<alphabet_t>;
        for (size_t i = 0; i < tmp.size(); ++i)
            seqan3::assign_rank_to(static_cast<rank_t>(i), tmp[i]);

        return tmp;
    }()};

    template <seqan3::semialphabet alphabet_t>
        requires (seqan3::alphabet_size<alphabet_t> <= 4)
    class snp_variant
    {
        uint32_t _value : 2;
        uint32_t _position : 30;

    public:
        snp_variant() = default;
        constexpr explicit snp_variant(uint32_t pos, alphabet_t value) noexcept :
            _value{seqan3::to_rank(value)},
            _position{pos}
        {
        }

        template <seqan3::cereal_output_archive output_archive_t>
        uint32_t save_minimal(output_archive_t const &) const
        {
            uint32_t tmp = _value;
            tmp = (tmp << 30) | _position; // make space for the position.
            return tmp;
        }

        template <seqan3::cereal_input_archive input_archive_t>
        void load_minimal(input_archive_t const &, uint32_t const &tmp)
        {
            _position = tmp & ((1 << 30) - 1);
            _value = (tmp >> 30);
        }

    private:
        constexpr friend uint32_t tag_invoke(std::tag_t<libjst::deletion>, snp_variant const &) noexcept
        {
            return 1;
        }

        constexpr friend std::span<alphabet_t const> tag_invoke(std::tag_t<libjst::insertion>, snp_variant const &me) noexcept
        {
            return {snp_value<alphabet_t>.data() + me._value, 1};
        }

        constexpr friend uint32_t tag_invoke(std::tag_t<libjst::position>, snp_variant const &me) noexcept
        {
            return me._position;
        }
    };
} // namespace libjst
