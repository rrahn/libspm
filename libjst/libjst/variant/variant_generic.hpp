// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an implementation for a generic sequence variant.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/concept.hpp>

#include <libjst/variant/concept.hpp>

namespace libjst
{
    template <seqan3::semialphabet alphabet_t>
    class generic_variant
    {
        std::vector<alphabet_t> _insertion{};
        uint32_t _position{};
        uint32_t _deletion{};

    public:
        constexpr generic_variant() = default;
        constexpr generic_variant(uint32_t const pos, std::vector<alphabet_t> ins, uint32_t const del) noexcept :
            _insertion{std::move(ins)},
            _position{pos},
            _deletion{del}
        {
        }

        template <seqan3::cereal_input_archive archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_insertion, _position, _deletion);
        }

        template <seqan3::cereal_output_archive archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_insertion, _position, _deletion);
        }

    private:
        constexpr friend uint32_t tag_invoke(std::tag_t<libjst::deletion>, generic_variant const &me) noexcept
        {
            return me._deletion;
        }

        constexpr friend std::span<alphabet_t const> tag_invoke(std::tag_t<libjst::insertion>, generic_variant const &me) noexcept
        {
            return {me._insertion};
        }

        constexpr friend uint32_t tag_invoke(std::tag_t<libjst::position>, generic_variant const &me) noexcept
        {
            return me._position;
        }
    };

} // namespace libjst
