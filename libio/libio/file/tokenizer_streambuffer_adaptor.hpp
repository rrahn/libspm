// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libio::tokenizer_streambuffer_adaptor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

namespace libio
{
    template <typename streambuffer_t>
    class tokenizer_streambuffer_adaptor : public std::ranges::view_base
    {
    public:
        using char_type = typename streambuffer_t::char_type;
        using traits_type = typename streambuffer_t::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;

    private:
        struct get_area_t : public streambuffer_t
        {
            using streambuffer_t::eback;
            using streambuffer_t::egptr;
            using streambuffer_t::gbump;
            using streambuffer_t::gptr;
            using streambuffer_t::setg;
            using streambuffer_t::uflow;
            using streambuffer_t::underflow;
        };

        class iterator;
        using sentinel = std::default_sentinel_t;

        get_area_t *_get_area{};

    public:
        tokenizer_streambuffer_adaptor() = default;

        constexpr explicit tokenizer_streambuffer_adaptor(streambuffer_t *stream_buffer) :
            _get_area{reinterpret_cast<get_area_t *>(stream_buffer)}
        {
            assert(_get_area != nullptr);
        }

        constexpr iterator begin() noexcept
        {
            assert(_get_area != nullptr);
            return iterator{this};
        }
        iterator begin() const noexcept = delete;

        constexpr sentinel end() noexcept
        {
            assert(_get_area != nullptr);
            return std::default_sentinel;
        }

        sentinel end() const noexcept = delete;
    };

    template <typename streambuffer_t>
    class tokenizer_streambuffer_adaptor<streambuffer_t>::iterator
    {
    private:
        using ger_area_iterator = char_type *;

        tokenizer_streambuffer_adaptor *_host{};
        int_type _eof_code{traits_type::eof()};
        ger_area_iterator _get_begin{nullptr};
        ger_area_iterator _get_end{nullptr};

    public:
        using value_type = std::span<char_type const>;
        using reference = std::span<char_type const>;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(tokenizer_streambuffer_adaptor *host) noexcept(noexcept(this->reset_get_area())) : _host{host}
        {
            reset_get_area();
        }

        constexpr reference operator*() const noexcept
        {
            return reference{_get_begin, _get_end};
        }

        constexpr iterator &operator++() noexcept(noexcept(this->reset_get_area()))
        {
            bump(_get_end - _get_begin);

            return *this;
        }

        constexpr void operator++(int) noexcept(noexcept(this->reset_get_area()))
        {
            ++(*this);
        }

        constexpr bool operator==([[maybe_unused]] std::default_sentinel_t const &rhs) const noexcept
        {
            return _eof_code == traits_type::eof();
        }

        constexpr void bump(difference_type offset) noexcept(noexcept(this->reset_get_area()))
        {
            _host->_get_area->gbump(offset);
            reset_get_area();
        }

    private:
        constexpr void reset_get_area() noexcept(noexcept(_host->_get_area->underflow()))
        {
            _eof_code = _host->_get_area->underflow(); // make sure enough data is available.
            _get_begin = _host->_get_area->gptr();
            _get_end = _host->_get_area->egptr();
        }
    };
} // namespace libio
