// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the delimitted token range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <iterator>
#include <iosfwd>
#include <span>

#include <libio/utility/tag_invoke.hpp>
#include <libio/utility/type_traits.hpp>

namespace libio
{
    namespace _stream_buffer {
        inline constexpr struct _cpo {

            template <typename tokenized_buffer_t>
                requires tag_invocable<_cpo, tokenized_buffer_t &>
            constexpr auto operator()(tokenized_buffer_t & buffer) const
                noexcept(is_nothrow_tag_invocable_v<_cpo, tokenized_buffer_t &>)
                -> tag_invoke_result_t<_cpo, tokenized_buffer_t &>
            {
                return libio::tag_invoke(_cpo{}, buffer);
            }

        } stream_buffer{};

    } // namespace _stream_buffer

    using _stream_buffer::stream_buffer;

    template <typename streambuffer_t>
    class tokenized_stream_buffer
    {
    public:
        using char_type = typename streambuffer_t::char_type;
        using traits_type = typename streambuffer_t::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;
    private:

        struct local_buffer_t : public streambuffer_t
        {
            using streambuffer_t::eback;
            using streambuffer_t::gptr;
            using streambuffer_t::egptr;
            using streambuffer_t::uflow;
            using streambuffer_t::underflow;
            using streambuffer_t::setg;
            using streambuffer_t::gbump;
        };

        using stop_fn_t = std::function<bool(char_type const)>;

        // multi-range iterator.
        class iterator;
        using sentinel = std::default_sentinel_t;

        local_buffer_t * _buffer{};
        stop_fn_t _stop_fn{};

    public:
        tokenized_stream_buffer() = default;

        template <typename delimiter_fn_t>
        constexpr tokenized_stream_buffer(streambuffer_t * stream_buffer, delimiter_fn_t delimiter_fn) :
            _buffer{reinterpret_cast<local_buffer_t *>(stream_buffer)},
            _stop_fn{std::move(delimiter_fn)}
        {
            assert(_buffer != nullptr);

            // on construction ensure that uflow is called and that eof is not set.
            int_type first = _buffer->uflow();
            if (first != traits_type::eof()) //
            {
                assert(_stop_fn(traits_type::to_char_type(first))); // expect to be at the begin delimiter of a token
            }
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

    private:

        template <typename this_t>
            requires std::same_as<std::remove_const_t<this_t>, tokenized_stream_buffer>
        friend constexpr auto tag_invoke(tag_t<libio::stream_buffer>, this_t & me) noexcept
            -> member_type_t<this_t &, local_buffer_t>
        {
            return *me._buffer;
        }

    };

    template <typename streambuffer_t>
    class tokenized_stream_buffer<streambuffer_t>::iterator
    {
    private:

        using buffer_iterator = char_type *;

        tokenized_stream_buffer * _host{};
        int_type _eof_char{traits_type::eof()};
        buffer_iterator _token_begin{nullptr};
        buffer_iterator _token_end{nullptr};

    public:

        // now we want to extend this to run on multiple ranges.
        using value_type = std::span<char_type const>;
        using reference = std::span<char_type const>;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(tokenized_stream_buffer * host) : _host{host}
        {
            update_stream_buffer();
        }

        constexpr reference operator*() const noexcept
        {
            return reference{_token_begin, _token_end};
        }

        constexpr iterator &operator++() noexcept(noexcept(libio::stream_buffer(*_host).sgetc()))
        {
            // advance buffer.
            libio::stream_buffer(*_host).gbump(_token_end - _token_begin);
            update_stream_buffer();

            return *this;
        }

        constexpr void operator++(int) noexcept(noexcept(libio::stream_buffer(*_host).sgetc()))
        {
            ++(*this);
        }

        constexpr bool operator==([[maybe_unused]] std::default_sentinel_t const &rhs) const noexcept
        {
            // either reload set the stream to eof, or token end points to something different then egptr()?
            return _eof_char == traits_type::eof();
        }
    private:

        constexpr void update_stream_buffer() noexcept(noexcept(libio::stream_buffer(*_host).sgetc())) {
            _eof_char = libio::stream_buffer(*_host).sgetc(); // only reloads the buffer if necessary, otherwise returns current symbol and does not reset the pointer.
            if (_eof_char != traits_type::eof()) {
                _token_begin = libio::stream_buffer(*_host).gptr();
                if (_token_begin == _token_end) { // not the same in the beginning.
                    _eof_char = traits_type::eof();
                } else {
                    _token_end = std::ranges::find_if(_token_begin, libio::stream_buffer(*_host).egptr(), _host->_stop_fn);
                }
            }
        }
    };
} // namespace libio
