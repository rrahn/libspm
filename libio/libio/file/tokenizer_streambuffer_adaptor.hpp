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
#include <span>

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
        tokenizer_streambuffer_adaptor(tokenizer_streambuffer_adaptor const &) = delete;
        tokenizer_streambuffer_adaptor(tokenizer_streambuffer_adaptor && other) noexcept :
            tokenizer_streambuffer_adaptor{}
        {
            using std::swap;
            swap(_get_area, other._get_area);
        }
        tokenizer_streambuffer_adaptor &operator=(tokenizer_streambuffer_adaptor const &) = delete;
        tokenizer_streambuffer_adaptor &operator=(tokenizer_streambuffer_adaptor && other) noexcept
        {
            _get_area = other._get_area;
            other._get_area = nullptr;
            return *this;
        }

        constexpr explicit tokenizer_streambuffer_adaptor(streambuffer_t *stream_buffer) :
            _get_area{reinterpret_cast<get_area_t *>(stream_buffer)}
        {
        }

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
        constexpr explicit iterator(tokenizer_streambuffer_adaptor *host) noexcept(noexcept(this->reset_get_area())) :
            _host{host}
        {
            if (_host->_get_area != nullptr)
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
            assert(offset <= _host->_get_area->in_avail());
            _host->_get_area->gbump(offset);
            reset_get_area();
        }

        constexpr void bump_with_restore(difference_type const offset) noexcept(noexcept(this->reset_get_area()))
        {
            assert(offset <= _host->_get_area->in_avail());
            auto gptr_curr = _host->_get_area->gptr();
            _host->_get_area->gbump(offset);
            if (_host->_get_area->gptr() == _get_end)  // reached end.
            {
                auto eback_original = _host->_get_area->eback();
                auto eback_new = std::ranges::copy(gptr_curr, _get_end, eback_original).out; // restore the skipped suffix in put back area
                _host->_get_area->setg(eback_new, _get_end, _get_end); // set the new dimension of the put back area.
                reset_get_area(); // reset the get area.
                assert(_get_begin == eback_new);  // make sure the underlying eback buffer was not removed.
                _get_begin = eback_original;
                _host->_get_area->setg(_get_begin, _get_begin, _get_end);
            }
            else // did not reach end.
            {
                reset_get_area();
            }
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
