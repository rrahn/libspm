// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the basic adapter for the seqan::Pattern objects.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <tuple>

#include <seqan/find.h>

#include <libcontrib/seqan/concept.hpp>
#include <libcontrib/seqan/container_adapter.hpp>

#include <libjst/matcher/concept.hpp>

namespace libjst
{
    template <typename derived_t>
    class seqan_pattern_base
    {
    private:

        friend derived_t;

        seqan_pattern_base() = default;

    public:

        // Note const is disabled since seqan use non-const pattern ;(
        template <std::ranges::viewable_range haystack_t, typename callback_t>
        constexpr void operator()(haystack_t && haystack, callback_t && callback) /*const*/ noexcept {
            using compatible_haystack_t = jst::contrib::seqan_container_t<std::views::all_t<haystack_t>>;

            if (std::ranges::empty(haystack)) {
                std::cout << "empty haystack\n";
                return;
            }

            compatible_haystack_t seqan_haystack =
                jst::contrib::make_seqan_container(std::views::all((haystack_t &&)haystack));

            seqan::Finder<compatible_haystack_t> finder{seqan_haystack};

            while (find_impl(finder, to_derived(this)->get_pattern())) {
                callback(finder);
            }
        }

    private:

        template <typename seqan_finder_t, typename seqan_pattern_t>
        constexpr bool find_impl(seqan_finder_t & finder, seqan_pattern_t && pattern) const noexcept {
            return std::apply([&] (auto && ...custom_args) {
                return seqan::find(finder, pattern, (decltype(custom_args)) custom_args...);
            }, to_derived(this)->custom_find_arguments());
        }

        constexpr auto & get_pattern() noexcept
        {
            return to_derived(this)->_pattern;
        }

        constexpr auto const & get_pattern() const noexcept
        {
            return to_derived(this)->_pattern;
        }

        constexpr auto custom_find_arguments() const noexcept {
            return std::tuple{};
        }

        static constexpr derived_t * to_derived(seqan_pattern_base * me) noexcept
        {
            return static_cast<derived_t *>(me);
        }

        static constexpr derived_t const * to_derived(seqan_pattern_base const * me) noexcept
        {
            return static_cast<derived_t const *>(me);
        }

        constexpr friend std::size_t tag_invoke(std::tag_t<libjst::window_size>, seqan_pattern_base const & me) noexcept {
            return seqan::length(seqan::needle(me.get_pattern()));
        }
    };

}  // namespace libjst
