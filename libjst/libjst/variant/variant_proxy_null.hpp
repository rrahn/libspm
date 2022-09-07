// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides nullified variant.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <std::integral position_t, typename coverage_t>
    class null_variant
    {
    private:
        position_t _position;
        coverage_t _coverage;
    public:

        template <std::integral _position_t, typename _coverage_t>
            requires (std::constructible_from<coverage_t, _coverage_t>)
        explicit null_variant(_position_t position, _coverage_t && coverage) :
            _position{position},
            _coverage{(_coverage_t &&)coverage}
        {}

    private:

        constexpr friend position_t tag_invoke(std::tag_t<libjst::position>, null_variant const & me) noexcept
        {
            return me._position;
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::insertion>, null_variant const &) noexcept
            -> libjst::variant_insertion_t<variant_t>
        {
            return libjst::variant_insertion_t<variant_t>{};
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::deletion>, null_variant const &) noexcept
            -> libjst::variant_deletion_t<variant_t>
        {
            return 0;
        }

        constexpr friend coverage_t const & tag_invoke(std::tag_t<libjst::coverage>, null_variant const &me) noexcept
        {
            return me._coverage;
        }
    };

    template <std::integral position_t, typename coverage_t>
    null_variant(position_t, coverage_t &&) -> null_variant<position_t, coverage_t>;
}  // namespace libjst
