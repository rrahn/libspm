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
#include <ranges>
#include <utility>

namespace libio
{
    // borrowed
    template <typename tokenizer_t>
    // requires tokenizer<tokenizer_t> // which is also a view
    class consume_tokenizer : public std::ranges::view_base
    {
    public:
        using char_type = typename std::remove_reference_t<tokenizer_t>::char_type;
        using traits_type = typename std::remove_reference_t<tokenizer_t>::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;

    private:
        using tokenizer_iterator = std::ranges::iterator_t<tokenizer_t>;
        class iterator;
        using sentinel = std::ranges::sentinel_t<tokenizer_t>;

        tokenizer_t _tokenizer{};
        tokenizer_iterator _cached_iter{};
        bool _called_begin{false};

    public:
        consume_tokenizer() = delete;
        constexpr explicit consume_tokenizer(tokenizer_t tokenizer) noexcept :
            _tokenizer{(tokenizer_t &&)tokenizer}
        {
        }
        consume_tokenizer(consume_tokenizer const &) = delete;
        consume_tokenizer(consume_tokenizer &&) = default;

        template <typename _tokenizer_t, typename ...args_t>
            requires std::constructible_from<_tokenizer_t, args_t...>
        constexpr explicit consume_tokenizer(std::in_place_type_t<_tokenizer_t>, args_t &&...args) noexcept :
            consume_tokenizer{_tokenizer_t{(args_t &&)args...}}
        {}

        ~consume_tokenizer()
        {
            begin(); // make sure begin on underlying iterator is called.
            while (_cached_iter != std::ranges::end(_tokenizer))
                ++_cached_iter;
        }

        constexpr iterator begin() noexcept
        {
            if (!_called_begin) // consumes but only once.
                _cached_iter = std::ranges::begin(_tokenizer);

            _called_begin = true;
            return iterator{this};
        }
        iterator begin() const noexcept = delete;

        constexpr sentinel end() noexcept
        {
            return std::ranges::end(_tokenizer);
        }

        sentinel end() const noexcept = delete;
    };

    template <typename tokenizer_t>
    consume_tokenizer(tokenizer_t &&) -> consume_tokenizer<tokenizer_t>;

    template <typename tokenizer_t, typename ...args_t>
    consume_tokenizer(std::in_place_type_t<tokenizer_t>, args_t &&...) -> consume_tokenizer<tokenizer_t>;

    template <typename tokenizer_t>
    class consume_tokenizer<tokenizer_t>::iterator
    {
    private:

        consume_tokenizer *_host{};

    public:
        using value_type = std::iter_value_t<tokenizer_iterator>; // return a new token abstraction
        using reference = std::iter_reference_t<tokenizer_iterator>;
        using difference_type = std::iter_difference_t<tokenizer_iterator>;
        using pointer = typename tokenizer_iterator::pointer;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(consume_tokenizer *host) noexcept : _host{host}
        {
        }

        constexpr reference operator*() const noexcept
        {
            return *(_host->_cached_iter);
        }

        constexpr iterator &operator++() noexcept
        {
            ++_host->_cached_iter;
            return *this;
        }

        constexpr void operator++(int) noexcept
        {
            ++(*this);
        }

        constexpr bool operator==(sentinel const &rhs) const noexcept
        {
            return _host->_cached_iter == rhs;
        }

        constexpr void bump(difference_type const offset) noexcept(noexcept(_host->_cached_iter.bump(offset)))
        {
            _host->_cached_iter.bump(offset);
        }
    };

} // namespace libio
