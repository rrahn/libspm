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

#include <libio/file/find_if.hpp>
#include <libio/file/until_tokenizer.hpp>

namespace libio
{
    template <typename tokenizer_t>
    class line_tokenizer : std::ranges::view_base //until_tokenizer<tokenizer_t>
    {
    private:
        using base_t = until_tokenizer<tokenizer_t>;

        // make more versatile by handling it as a multicharacter matcher.
        static constexpr auto is_newline = [] (char const c) { return c == '\n' || c == '\r'; }; //seqan3::is_char<'\n'> || seqan3::is_char<'\r'>;

        class iterator;
        using sentinel = std::ranges::sentinel_t<tokenizer_t>;

        // may create _until tokenizer and store original tokenizer
        tokenizer_t _tokenizer;
        // until_tokenizer<tokenizer_t> _until_tokenizer;

    public:
        using char_type = typename std::remove_reference_t<tokenizer_t>::char_type;
        using traits_type = typename std::remove_reference_t<tokenizer_t>::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;

        line_tokenizer() = delete;
        constexpr explicit line_tokenizer(tokenizer_t tokenizer) noexcept : //base_t{(tokenizer_t &&)tokenizer, is_newline}
            _tokenizer{(tokenizer_t &&)tokenizer}//,
            // _until_tokenizer{_tokenizer, is_newline}
        {
        }
        line_tokenizer(line_tokenizer const &) = delete;
        line_tokenizer(line_tokenizer &&) = default;

        constexpr iterator begin() noexcept
        {
            return iterator{this};
        }
        void begin() const noexcept = delete;

        constexpr sentinel end() noexcept
        {
            return std::ranges::end(_tokenizer);
        }

        void end() const noexcept = delete;
    };

    template <typename tokenizer_t>
    line_tokenizer(tokenizer_t &&) -> line_tokenizer<tokenizer_t>;

   // Closure object
    struct line_token
    {
        template <typename tokenizer_t>
        constexpr auto operator()(tokenizer_t && tokenizer) const noexcept
        {
            return line_tokenizer{(tokenizer_t &&) tokenizer};
        }
    };

    template <typename tokenizer_t>
    class line_tokenizer<tokenizer_t>::iterator
    {
    private:
        using tokenizer_iterator = std::ranges::iterator_t<tokenizer_t>;
        using get_area_t = std::iter_reference_t<tokenizer_iterator>;
        using get_area_iterator = std::ranges::iterator_t<get_area_t>;

        line_tokenizer *_host{};
        get_area_iterator _get_begin{};
        get_area_iterator _get_end{};
        get_area_iterator _token_end{};
        tokenizer_iterator _it{};
        bool _found_token_end{false};

    public:
        using value_type = get_area_t; // return a new token abstraction
        using reference = get_area_t;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(line_tokenizer *host) noexcept : _host{host}
        {
            // this is the area until the end of the tokenizer
            _it = std::ranges::begin(_host->_tokenizer); // initialise underlying stream buffer
            if (_it != std::ranges::end(_host->_tokenizer))
            {
                underflow();
            }
                // reset_get_area();
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
            return _it == rhs || _get_begin == _get_end;
        }

        constexpr void bump(difference_type const offset) noexcept(noexcept(_it.bump(offset)))
        {
            assert(offset <= _get_end - _get_begin);
            _it.bump(offset);
            reset_get_area(offset);
        }

        constexpr void bump_with_restore(difference_type const offset) noexcept(noexcept(_it.bump_with_restore(offset)))
        {
            _it.bump_with_restore(offset);
            reset_get_area(offset);
        }

    private:
        constexpr void reset_get_area(difference_type const offset) noexcept
        {
            _get_begin += offset;
            if (_get_begin == _get_end)  // set to end condition.
            {
                if (_found_token_end) // found token inside current get area.
                {
                    _it.bump(_token_end - _get_end); // consume new line characters.
                    _get_begin += _token_end - _get_end;
                    _get_end += _token_end - _get_end;
                }
                else
                {
                    underflow(); // underflow and call again.
                    // case 5: begin == end && token end found => we need to signal that we are not at end
                    if (_get_begin == _get_end && _found_token_end)  // set to end condition.
                    {
                        _it.bump(_token_end - _get_end); // consume new line characters.
                        _get_begin += _token_end - _get_end;
                        _get_end += _token_end - _get_end;
                    }
                    // case 1: eof stream -> next time we compare we reached eof.
                    // case 2: begin < end && token end not found -> we can just continue searching
                    // case 3: begin < end && token end found -> we can just continue searching
                    // case 4: begin == end  && token end not found => same as case 1: -> we stop
                }
            }
        }

        constexpr void underflow() noexcept
        {
            get_area_t const &get_area = *_it;
            _get_begin = get_area.begin(); // get new area.
            _get_end = std::ranges::find_if(get_area, is_newline);
            _token_end = std::ranges::find_if_not(_get_end, get_area.end(), is_newline);
            _found_token_end = _get_end != _token_end;
        }
    };
} // namespace libio
