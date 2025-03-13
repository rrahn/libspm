// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides restorable adaption of myers matcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libspm/matcher/seqan_pattern_base.hpp>
#include <libspm/matcher/seqan_restorable_pattern.hpp>

namespace seqan2 {

    template <typename needle_t>
    class Pattern<needle_t, spm::Restorable<MyersUkkonenGlobal>> :
        public Pattern<needle_t, MyersUkkonenGlobal>
    {
    private:
        using base_type = Pattern<needle_t, MyersUkkonenGlobal>;

        bool _first_find{true};

    public:

        using state_type = typename base_type::TPatternState;

        Pattern() = default;

        template <std::ranges::viewable_range _needle_t, std::unsigned_integral error_count_t>
            requires (!std::same_as<_needle_t, Pattern>)
        constexpr explicit Pattern(_needle_t && needle, error_count_t const max_error_count = 0u) : base_type{}
        {
            if (!std::ranges::empty(needle)) {
                setHost(to_base(), (_needle_t &&)needle);
                setScoreLimit(to_base(), -static_cast<int32_t>(max_error_count));
                _patternFirstInit(to_base(), seqan2::needle(to_base()));
                _patternInit(to_base(), to_state(), std::ignore);
            }
        }

        template <typename finder_t>
        constexpr bool operator()(finder_t & finder) noexcept {
            using haystack_t = typename Haystack<finder_t>::Type;
            using haystack_size_t = typename Size<haystack_t>::Type;

            if (empty(to_base().data_host) || !initialise(finder)) return false;

            haystack_size_t haystack_length = std::min<haystack_size_t>(length(haystack(finder)),
                                                                        this->needleSize - scoreLimit(to_state()) + 1);

            if (is_short())
                return _findMyersSmallPatterns(finder, to_base(), to_state(), haystack_length);
            else
                return _findMyersLargePatterns(finder, to_base(), to_state(), haystack_length);
        }

        constexpr state_type const & capture() const noexcept {
            return to_state();
        }

        constexpr void restore(state_type state) noexcept {
            to_state() = std::move(state);
        }

    protected:
        constexpr bool is_short() const noexcept {
            return this->largePattern == nullptr;
        }

    private:

        template <typename finder_t>
        constexpr bool initialise(finder_t & finder) noexcept {
            if (!_first_find && !empty(finder)) {
                if (atEnd(finder))
                    return false;
                goNext(finder);
            }
            _first_find = false;
            _finderSetNonEmpty(finder);
            return true;
        }

        constexpr base_type & to_base() noexcept {
            return static_cast<base_type &>(*this);
        }

        constexpr base_type const & to_base() const noexcept {
            return static_cast<base_type const &>(*this);
        }

        constexpr state_type & to_state() noexcept {
            return static_cast<state_type &>(*this);
        }

        constexpr state_type const & to_state() const noexcept {
            return static_cast<state_type const &>(*this);
        }
    };

    template <typename finder_t, typename needle_t>
    constexpr bool find(finder_t & finder, Pattern<needle_t, spm::Restorable<MyersUkkonenGlobal>> & pattern) {
        return pattern(finder);
    }

} // namespace seqan2

namespace spm
{

    template <std::ranges::random_access_range needle_t>
    class restorable_myers_prefix_matcher : public seqan_pattern_base<restorable_myers_prefix_matcher<needle_t>>
    {
    private:

        using base_t = seqan_pattern_base<restorable_myers_prefix_matcher<needle_t>>;

        friend base_t;

        using compatible_needle_type = spm::seqan_container_t<needle_t>;
        using pattern_type = seqan2::Pattern<compatible_needle_type, Restorable<seqan2::MyersUkkonenGlobal>>;

        pattern_type _pattern{};

    public:

        using state_type = typename pattern_type::state_type;

        restorable_myers_prefix_matcher() = delete;
        template <std::ranges::viewable_range _needle_t, std::unsigned_integral error_count_t>
            requires (!std::same_as<_needle_t, restorable_myers_prefix_matcher> &&
                       std::constructible_from<compatible_needle_type, _needle_t>)
        explicit restorable_myers_prefix_matcher(_needle_t && needle, error_count_t const error_count) :
            _pattern{spm::make_seqan_container(std::views::all((_needle_t &&) needle)), error_count}
        {}

        constexpr state_type const & capture() const noexcept {
            return _pattern.capture();
        }

        constexpr void restore(state_type state) noexcept {
            _pattern.restore(std::move(state));
        }

    private:

        constexpr pattern_type & get_pattern() noexcept {
            return _pattern;
        }

        constexpr friend std::size_t tag_invoke(std::tag_t<window_size>, restorable_myers_prefix_matcher const & me) noexcept {
            return spm::window_size(static_cast<base_t const &>(me)) - seqan2::scoreLimit(me._pattern.capture());
        }
    };

    template <std::ranges::viewable_range needle_t, std::unsigned_integral error_count_t>
    restorable_myers_prefix_matcher(needle_t &&, error_count_t) -> restorable_myers_prefix_matcher<std::views::all_t<needle_t>>;

}  // namespace spm
