// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an adapter to make the shiftor online pattern matching algorithm work with the JST.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libspm/matcher/seqan_pattern_base.hpp>

namespace spm
{

    template <std::ranges::random_access_range needle_t>
    class shiftor_matcher : public seqan_pattern_base<shiftor_matcher<needle_t>>
    {
    private:

        friend seqan_pattern_base<shiftor_matcher<needle_t>>;

        using compatible_needle_type = spm::seqan_container_t<needle_t>;
        using pattern_type = seqan2::Pattern<compatible_needle_type, seqan2::ShiftOr>;

        pattern_type _pattern{};

    public:

        shiftor_matcher() = delete;
        template <std::ranges::viewable_range _needle_t>
            requires (!std::same_as<_needle_t, shiftor_matcher> &&
                       std::constructible_from<compatible_needle_type, _needle_t>)
        explicit shiftor_matcher(_needle_t && needle) :
            _pattern{spm::make_seqan_container(std::views::all((_needle_t &&) needle))}
        {}
    };

    template <std::ranges::viewable_range needle_t>
    shiftor_matcher(needle_t &&) -> shiftor_matcher<std::views::all_t<needle_t>>;

}  // namespace spm
