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

#include <ranges>

#include <seqan3/utility/char_operations/predicate.hpp>

#include <libio/file/until_tokenizer.hpp>

namespace libio
{
    template <typename tokenizer_t>
    class line_tokenizer : public until_tokenizer<tokenizer_t>
    {
    private:
        using base_t = until_tokenizer<tokenizer_t>;

        // make more versatile by handling it as a multicharacter matcher.
        static constexpr auto is_newline = seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;
    public:
        using typename base_t::char_type;
        using typename base_t::traits_type;
        using typename base_t::int_type;
        using typename base_t::pos_type;
        using typename base_t::off_type;

        line_tokenizer() = default;

        constexpr explicit line_tokenizer(tokenizer_t tokenizer) noexcept : base_t{std::move(tokenizer), is_newline}
        {
        }

        // ~line_buffer() {
        //     // consume remaining part
        //     while ()
        // }

        constexpr std::ranges::iterator_t<base_t> begin() noexcept
        {
            return base_t::begin();
        }
        void begin() const noexcept = delete;

        constexpr std::ranges::sentinel_t<base_t> end() noexcept
        {
            return base_t::end();
        }

        void end() const noexcept = delete;
    };

    // template <typename tokenizer_t>
    // class line_tokenizer<tokenizer_t>::iterator
    // {
    // private:
    //     static constexpr auto is_newline = seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;

    //     using char_type = typename tokenizer_t::char_type;
    //     using traits_type = typename tokenizer_t::traits_type;

    //     using buffer_iterator = std::ranges::iterator_t<tokenizer_t>;
    //     using buffer_span_t = std::iter_reference_t<buffer_iterator>;
    //     using buffer_span_iterator = std::ranges::iterator_t<buffer_span_t>;

    //     line_tokenizer *_host{};
    //     buffer_span_iterator _segment_begin{};
    //     buffer_span_iterator _segment_end{};
    //     buffer_iterator _it{};
    //     bool _found_newline{false};

    // public:
    //     using value_type = buffer_span_t;
    //     using reference = buffer_span_t;
    //     using difference_type = std::ptrdiff_t;
    //     using pointer = void;
    //     using category = std::input_iterator_tag;

    //     iterator() = default;
    //     constexpr explicit iterator(line_tokenizer *host) noexcept : _host{host}
    //     {
    //         _it = std::ranges::begin(_host->_buffer); // defines the range of the iterator.
    //         _segment_end = (*_it).begin();
    //         set_buffer();
    //     }

    //     constexpr reference operator*() const noexcept
    //     {
    //         return reference{_segment_begin, _segment_end};
    //     }

    //     constexpr iterator &operator++() noexcept
    //     {
    //         set_buffer();
    //         return *this;
    //     }

    //     constexpr void operator++(int) noexcept
    //     {
    //         ++(*this);
    //     }

    //     constexpr bool operator==(sentinel const &rhs) const noexcept
    //     {
    //         return _it == rhs || _found_newline;
    //     }

    // private:
    //     constexpr void set_buffer() noexcept(noexcept(++_it))
    //     {
    //         _segment_begin = _segment_end; // reset
    //         _segment_end = std::ranges::find_if(_segment_begin, (*_it).end(), is_newline);

    //         if (_segment_begin == _segment_end)
    //         {
    //             if (_segment_begin == (*_it).end())
    //             {
    //                 ++_it; // reload underlying buffer.
    //                 _segment_end = (*_it).begin();
    //             }
    //             else
    //             { // reset the buffer: better if we can do this here.
    //                 // we must let the underlying iterator do this.
    //                 // we need some version to manifest that we gbumped?
    //                 // how many fields did we go?
    //                 libio::stream_buffer(*_host).gbump(std::ranges::distance((*_it).begin(), _segment_end));
    //                 // libio::stream_buffer(*_host).setg(libio::stream_buffer(*_host).eback(),
    //                 //                                   _segment_begin,
    //                 //                                   libio::stream_buffer(*_host).egptr());
    //                 _found_newline = true;
    //             }
    //         }

    //         // // we should find the end?
    //         // _token_begin =

    //         // _eof_char = libio::stream_buffer(*_host).sgetc(); // only reloads the buffer if necessary, otherwise returns current symbol and does not reset the pointer.
    //         // if (_eof_char != traits_type::eof()) {
    //         //     _token_begin = libio::stream_buffer(*_host).gptr();
    //         //     if (_token_begin == _token_end) { // not the same in the beginning.
    //         //         bool is_cr = seqan3::is_char<'\r'>(*_token_begin); // read current element
    //         //         bool is_lf = seqan3::is_char<'\n'>(traits_type::to_char_type(libio::stream_buffer(*_host).snextc())); // maybe call uflow if at end of symbol.
    //         //         if (is_cr && is_lf) libio::stream_buffer(*_host).sbumpc();
    //         //         _eof_char = traits_type::eof(); // set end of stream
    //         //     } else {
    //         //         _token_end = std::ranges::find_if(_token_begin, libio::stream_buffer(*_host).egptr(), is_newline);
    //         //     }
    //         // }
    //     }
    // };

} // namespace libio
