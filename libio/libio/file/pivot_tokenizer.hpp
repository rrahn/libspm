// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides pivot tokenizer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <algorithm>
#include <functional>
#include <ranges>
#include <string_view>
#include <span>

namespace libio
{

    template <typename char_t, size_t N>
    class pivot_matcher
    {
    public:
        using needle_type = std::span<char_t const, N>;
    private:

        needle_type _needle{std::array<char_t, N>{}};

    public:
        pivot_matcher() = default;

        constexpr explicit pivot_matcher(char_t const(&needle)[N + 1]) noexcept : _needle{&needle[0], N}
        {
        }

        needle_type const & needle() const noexcept
        {
            return _needle;
        }

        constexpr std::span<const char_t> operator()(std::span<const char_t> const & haystack) const noexcept
        {
            using iter_t = std::ranges::iterator_t<std::span<const char_t>>;
            iter_t _hit_end = std::ranges::begin(haystack);
            size_t match_count{};
            while (_hit_end != std::ranges::end(haystack) && match_count < std::ranges::size(needle()))
            {
                if (*_hit_end != needle()[match_count])
                {
                    match_count = 0; // reset to 0.
                    _hit_end = std::ranges::find(_hit_end, std::ranges::end(haystack), _needle[0]);
                    continue;
                }
                ++match_count;
                ++_hit_end;
            }
            assert(static_cast<std::ptrdiff_t>(match_count) <= (_hit_end - std::ranges::begin(haystack)));
            return std::span{_hit_end - match_count, match_count};
        }
    };

    template <typename char_t, size_t N>
    pivot_matcher(char_t const(&)[N]) -> pivot_matcher<char_t, N - 1>;

    // borrowed
    template <typename tokenizer_t, typename matcher_t>
    // requires tokenizer<tokenizer_t> // which is also a view
    class pivot_tokenizer : public std::ranges::view_base
    {
    public:
        using char_type = typename std::remove_reference_t<tokenizer_t>::char_type;
        using traits_type = typename std::remove_reference_t<tokenizer_t>::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;

    private:

        class iterator;
        using sentinel = std::ranges::sentinel_t<tokenizer_t>;

        tokenizer_t _tokenizer;
        matcher_t _matcher{};

    public:
        pivot_tokenizer() = delete;

        constexpr explicit pivot_tokenizer(tokenizer_t tokenizer, matcher_t matcher) noexcept :
            _tokenizer{(tokenizer_t &&)tokenizer},
            _matcher{std::move(matcher)}
        {
        }
        pivot_tokenizer(pivot_tokenizer const &) = delete;
        pivot_tokenizer(pivot_tokenizer &&) = default;

        constexpr iterator begin() noexcept
        {
            return iterator{this};
        }
        iterator begin() const noexcept = delete;

        constexpr sentinel end() noexcept
        {
            return std::default_sentinel;
        }

        sentinel end() const noexcept = delete;
    };

    template <typename tokenizer_t, typename matcher_t>
    pivot_tokenizer(tokenizer_t &&, matcher_t const &) -> pivot_tokenizer<tokenizer_t, matcher_t>;

    // Closure object
    template <typename matcher_t>
    struct pivot_token
    {
        matcher_t _matcher{}; // the matcher is now const

        template <typename tokenizer_t>
        constexpr auto operator()(tokenizer_t && tokenizer) const noexcept
        {
            // using pivot_tokenizer_t = pivot_tokenizer<tokenizer_t, pivot_matcher_t>;
            return pivot_tokenizer{(tokenizer_t &&) tokenizer, std::move(_matcher)};
        }
    };

    template <typename matcher_t>
    pivot_token(matcher_t const &) -> pivot_token<matcher_t>;

    template <typename tokenizer_t, typename matcher_t>
    class pivot_tokenizer<tokenizer_t, matcher_t>::iterator
    {
    private:
        using tokenizer_iterator = std::ranges::iterator_t<tokenizer_t>;
        using get_area_t = std::iter_reference_t<tokenizer_iterator>;
        using get_area_iterator = std::ranges::iterator_t<get_area_t>;

        pivot_tokenizer *_host{};
        get_area_iterator _get_begin{};
        get_area_t _hit_area{};
        tokenizer_iterator _it{};

    public:
        using value_type = get_area_t; // return a new token abstraction
        using reference = get_area_t;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(pivot_tokenizer *host) noexcept : _host{host}
        {
            _it = std::ranges::begin(_host->_tokenizer); // initialise underlying stream buffer
            if (_it != std::ranges::end(_host->_tokenizer))
            {
                auto const & get_area = *_it;
                _get_begin = get_area.begin();
                _hit_area = _host->_matcher(get_area); // searching for the hit
            }
        }

        constexpr reference operator*() const noexcept
        {
            return reference{_get_begin, _hit_area.begin()};
        }

        constexpr iterator &operator++() noexcept
        {
            assert(_get_begin <= _hit_area.begin());
            bump(_hit_area.begin() - _get_begin);

            return *this;
        }

        constexpr void operator++(int) noexcept
        {
            ++(*this);
        }

        constexpr bool operator==(sentinel const &rhs) const noexcept
        {
            return _it == rhs || _get_begin == _hit_area.begin();
        }

        constexpr void bump(difference_type const offset) noexcept(noexcept(_it.bump(offset)))
        {
            assert(offset <= _hit_area.begin() - _get_begin); // offset must be smaller.
            _it.bump(offset); // update all underlying iterators to point to current range.
            // this might also set the underlying buffer to the end.
            if (_it != std::ranges::end(_host->_tokenizer))
            {   // case 1: hit_area can only be empty if it was the end of the get area.
                if (std::ranges::empty(_hit_area)) // no hit found in last get_area, must have been updated already
                {
                    auto const &get_area = *_it;
                    _get_begin = get_area.begin();
                    _hit_area = _host->_matcher(get_area);
                } // case 2: found full hit => implies that bump was not calling underflow because buffer not saturated
                else if (std::ranges::size(_hit_area) == std::ranges::size(_host->_matcher.needle()))
                {
                    _get_begin += offset; // also update current _get_begin
                    assert(_get_begin <= _hit_area.begin()); // also signal end
                } // case 3: found partial hit at end => implies bump did not caused call to underfow
                else
                {
                    assert(std::ranges::size(_hit_area) <= std::ranges::size(_host->_matcher.needle()));
                    assert(_hit_area.begin() == (*_it).begin()); // underlying is set to begin of hit.
                    assert(_hit_area.end() == (*_it).end()); // bump did not overflow yet.
                    bump_with_restore(std::ranges::ssize(_hit_area));
                }
            }
        }

        constexpr void bump_with_restore(difference_type const offset) noexcept(noexcept(_it.bump_with_restore(offset)))
        {
            _it.bump_with_restore(offset); // this delegates the update to the underlying buffer.
            get_area_t const &get_area = *_it;
            _get_begin = get_area.begin();
            _hit_area = _host->_matcher(get_area);
        }
    };

} // namespace libio
