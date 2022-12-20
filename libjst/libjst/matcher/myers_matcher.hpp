// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an adapter to make the myers online pattern matching algorithm work with the JST.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/matcher/seqan_pattern_base.hpp>

namespace libjst
{
    template <std::ranges::random_access_range needle_t>
    class myers_matcher : public seqan_pattern_base<myers_matcher<needle_t>>
    {
    private:

        using base_t = seqan_pattern_base<myers_matcher<needle_t>>;

        friend base_t;

        using compatible_needle_type = jst::contrib::seqan_container_t<needle_t>;
        using pattern_type = seqan::Pattern<compatible_needle_type, seqan::Myers<>>;

        pattern_type _pattern{};
        int32_t _min_score{};

    public:

        myers_matcher() = delete;
        template <std::ranges::viewable_range _needle_t>
            requires (!std::same_as<_needle_t, myers_matcher> &&
                       std::constructible_from<compatible_needle_type, _needle_t>)
        explicit myers_matcher(_needle_t && needle, std::size_t max_error_count = 0) :
            _pattern{jst::contrib::make_seqan_container(std::views::all((_needle_t &&) needle))},
            _min_score{-static_cast<int32_t>(max_error_count)}
        {}

    private:

        constexpr auto custom_find_arguments() const noexcept {
            return std::tuple{_min_score};
        }

        constexpr friend std::size_t tag_invoke(std::tag_t<window_size>, myers_matcher const & me) noexcept {
            return libjst::window_size(static_cast<base_t const &>(me)) - me._min_score;
        }
    };

    template <std::ranges::viewable_range needle_t>
    myers_matcher(needle_t &&) -> myers_matcher<std::views::all_t<needle_t>>;

    template <std::ranges::viewable_range needle_t>
    myers_matcher(needle_t &&, std::size_t) -> myers_matcher<std::views::all_t<needle_t>>;

}  // namespace libjst
