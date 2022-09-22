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

// #include <libio/file/tokenized_stream_buffer.hpp>

namespace libio
{
    template <typename tokenizer_t>
    class segment_tokenizer : public std::ranges::view_base
    {
    public:
        using char_type = typename std::remove_reference_t<tokenizer_t>::char_type;
        using traits_type = typename std::remove_reference_t<tokenizer_t>::traits_type;
        using int_type = typename traits_type::int_type;
        using pos_type = typename traits_type::pos_type;
        using off_type = typename traits_type::off_type;
    private:
        // faster options to tokenize a line?

        class iterator;
        using sentinel = std::ranges::sentinel_t<tokenizer_t>;
        using fn_t = std::function<bool(char_type const)>;

        tokenizer_t _tokenizer;
        fn_t _segment_fn{};

    public:
        segment_tokenizer() = delete;
        template <typename segment_fn_t>
        constexpr explicit segment_tokenizer(tokenizer_t tokenizer, segment_fn_t segment_fn) noexcept
            :  _tokenizer{(tokenizer_t &&)tokenizer}, _segment_fn{std::move(segment_fn)}
        {
        }
        segment_tokenizer(segment_tokenizer const &) = delete;
        segment_tokenizer(segment_tokenizer &&) = default;

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

    template <typename tokenizer_t, typename segment_fn_t>
    segment_tokenizer(tokenizer_t &&, segment_fn_t &&) -> segment_tokenizer<tokenizer_t>;

    template <typename tokenizer_t>
    class segment_tokenizer<tokenizer_t>::iterator
    {
    private:

        using tokenizer_iterator = std::ranges::iterator_t<tokenizer_t>;
        using get_area_t = std::iter_reference_t<tokenizer_iterator>;
        using get_area_iterator = std::ranges::iterator_t<get_area_t>;

        segment_tokenizer *_host{};
        get_area_iterator _get_begin{};
        get_area_iterator _get_next{};
        get_area_iterator _get_end{};
        tokenizer_iterator _it{};

    public:

        using value_type = get_area_t;
        using reference = get_area_t;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using category = std::input_iterator_tag;

        iterator() = default;
        constexpr explicit iterator(segment_tokenizer * host) noexcept : _host{host}
        {
            _it = std::ranges::begin(_host->_tokenizer); // defines the range of the iterator.
            if (_it != std::ranges::end(_host->_tokenizer))
                reset_get_area();
        }

        constexpr reference operator*() const noexcept
        {
            return reference{_get_next, _get_end}; // maybe empty
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

        constexpr bool operator==(sentinel const & rhs) const noexcept
        {
            return _it == rhs;
        }

        constexpr void bump(difference_type const offset) noexcept(noexcept(_it.bump(1)))
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

        constexpr void reset_get_area() noexcept(noexcept(++_it)) {
            // underlying it pointer has not changed!
            get_area_t const & get_area = *_it;
            _get_begin = get_area.begin();
            _get_next = std::ranges::find_if(_get_begin, get_area.end(), _host->_segment_fn);
            _get_end = std::ranges::find_if_not(_get_next, get_area.end(), _host->_segment_fn);
        }
    };

} // namespace libio
