// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides restorable adaption of shiftor matcher.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/matcher/seqan_pattern_base.hpp>
#include <libjst/matcher/seqan_restorable_pattern.hpp>

namespace seqan {
    template <typename needle_t>
    class Pattern<needle_t, libjst::Restorable<ShiftOr>> : public Pattern<needle_t, ShiftOr>
    {
    private:
        using base_type = Pattern<needle_t, ShiftOr>;

        bool _first_find{true};

    public:

        using state_type = decltype(base_type::prefSufMatch);

        Pattern() = default;
        constexpr explicit Pattern(needle_t const & needle) : base_type{needle}
        {}

        template <typename finder_t>
        constexpr bool operator()(finder_t & finder) noexcept {
            initialise(finder);
            if (is_short())
                return seqan::_findShiftOrSmallNeedle(finder, to_base(*this));
            else
                return seqan::_findShiftOrLargeNeedle(finder, to_base(*this));
        }

        constexpr state_type const & capture() const noexcept {
            return this->prefSufMatch;
        }

        constexpr void restore(state_type state) noexcept {
            this->prefSufMatch = std::move(state);
        }

    protected:

        template <typename finder_t>
        constexpr void initialise(finder_t & finder) noexcept {
            if (_first_find) {
                _patternInit(to_base(*this));
                _first_find = false;
            }
            _setFinderLength(finder, this->needleLength);
            _finderSetNonEmpty(finder);
            setPosition(finder, endPosition(finder));
        }

        constexpr bool is_short() const noexcept {
            return this->blockCount == 1;
        }

        private:

            static constexpr base_type & to_base(Pattern & me) noexcept {
                return static_cast<base_type &>(me);
            }

            static constexpr base_type const & to_base(Pattern const & me) noexcept {
                return static_cast<base_type const &>(me);
            }
    };

    template <typename finder_t, typename needle_t>
    constexpr bool find(finder_t & finder, Pattern<needle_t, libjst::Restorable<ShiftOr>> & pattern) {
        return pattern(finder);
    }
} // namespace seqan

namespace libjst
{

    template <std::ranges::random_access_range needle_t>
    class restorable_shiftor_matcher : public seqan_pattern_base<restorable_shiftor_matcher<needle_t>>
    {
    private:

        friend seqan_pattern_base<restorable_shiftor_matcher<needle_t>>;

        using compatible_needle_type = jst::contrib::seqan_container_t<needle_t>;
        using pattern_type = seqan::Pattern<compatible_needle_type, Restorable<seqan::ShiftOr>>;

        pattern_type _pattern{};

    public:

        using state_type = typename pattern_type::state_type;

        restorable_shiftor_matcher() = delete;
        template <std::ranges::viewable_range _needle_t>
            requires (!std::same_as<_needle_t, restorable_shiftor_matcher> &&
                       std::constructible_from<compatible_needle_type, _needle_t>)
        explicit restorable_shiftor_matcher(_needle_t && needle) :
            _pattern{jst::contrib::make_seqan_container(std::views::all((_needle_t &&) needle))}
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
    };

    template <std::ranges::viewable_range needle_t>
    restorable_shiftor_matcher(needle_t &&) -> restorable_shiftor_matcher<std::views::all_t<needle_t>>;

}  // namespace libjst
