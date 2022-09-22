// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implements the line buffer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <algorithm>
#include <functional>
#include <ranges>

namespace libio
{

    template <typename tokenizer_t>
    // requires tokenizer<tokenizer_t> // which is also a view
    class until_tokenizer : public std::ranges::view_base
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
        std::function<bool(char_type const)> _until_fn{};

    public:
        until_tokenizer() = delete;
        template <typename until_token_t>
        constexpr explicit until_tokenizer(tokenizer_t tokenizer, until_token_t &&until_fn) noexcept :
            _tokenizer{(tokenizer_t &&)tokenizer},
            _until_fn{until_fn}
        {
        }
        until_tokenizer(until_tokenizer const &) = delete;
        until_tokenizer(until_tokenizer &&) = default;

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

    template <typename tokenizer_t, typename until_token_t>
    until_tokenizer(tokenizer_t &&, until_token_t &&) -> until_tokenizer<tokenizer_t>;

    // Closure object
    template <typename predicate_t>
    struct until_token
    {
        predicate_t _predicate{};

        template <typename tokenizer_t>
        constexpr auto operator()(tokenizer_t && tokenizer) const noexcept
        {
            return until_tokenizer{(tokenizer_t &&) tokenizer, std::move(_predicate)};
        }
    };

    template <typename tokenizer_t>
    class until_tokenizer<tokenizer_t>::iterator
    {
    private:
        using tokenizer_iterator = std::ranges::iterator_t<tokenizer_t>;
        using get_area_t = std::iter_reference_t<tokenizer_iterator>;
        using get_area_iterator = std::ranges::iterator_t<get_area_t>;

        until_tokenizer *_host{};
        get_area_iterator _get_begin{};
        get_area_iterator _get_end{};
        tokenizer_iterator _it{};

    public:
        using value_type = get_area_t; // return a new token abstraction
        using reference = get_area_t;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(until_tokenizer *host) noexcept : _host{host}
        {
            _it = std::ranges::begin(_host->_tokenizer); // initialise underlying stream buffer

            while (_it != std::ranges::end(_host->_tokenizer))
            {
                auto get_area = *_it;
                if (auto initial = std::ranges::find_if_not(get_area, _host->_until_fn); initial != get_area.end())
                {
                    bump(initial - get_area.begin()); // found pattern so we set to that symbol.
                    break;
                }
                ++_it; // bumps entire lower range.
            }
        }

        constexpr reference operator*() const noexcept
        {
            return reference{_get_begin, _get_end};
        }

        constexpr iterator &operator++() noexcept
        {
            bump(_get_end - _get_begin);
            return *this;
        }

        constexpr void operator++(int) noexcept
        {
            ++(*this);
        }

        constexpr bool operator==(sentinel const &rhs) const noexcept
        {
            return _it == rhs || _host->_until_fn(*_get_begin);
        }

        constexpr void bump(difference_type const offset) noexcept(noexcept(_it.bump(offset)))
        {
            _it.bump(offset);
            reset_get_area();
        }

        constexpr void bump_with_restore(difference_type const offset) noexcept(noexcept(_it.bump_with_restore(offset)))
        {
            _it.bump_with_restore(offset);
            reset_get_area();
        }

    private:
        constexpr void reset_get_area() noexcept
        {
            get_area_t const &get_area = *_it;
            _get_begin = get_area.begin();
            _get_end = std::ranges::find_if(get_area, _host->_until_fn);
        }
    };

} // namespace libio
