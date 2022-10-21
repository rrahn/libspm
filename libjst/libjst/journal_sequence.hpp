// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal_sequence.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

namespace libjst
{

// TODO: rename to journaled_sequence with sequence interface.
 // uses internally the journal, which can be specified with some other journal type if we want to.
template <typename journal_t>
class journal_sequence
{
private:

    journal_t const * _journal{};

    class iterator;
public:

    journal_sequence() = default;
    explicit journal_sequence(journal_t const & journal) : _journal{std::addressof(journal)}
    {}

    size_t size() const noexcept
    {
        assert(_journal != nullptr);
        return _journal->sequence_size();
    }

    bool empty() const noexcept
    {
        return size() == 0;
    }

    iterator begin() const & noexcept
    {
        assert(_journal != nullptr);
        return iterator{*_journal, 0};
    }

    iterator begin() && noexcept
    {
        assert(_journal != nullptr);
        return iterator{*_journal, 0};
    }

    iterator end() const & noexcept
    {
        assert(_journal != nullptr);
        return iterator{*_journal, size()};
    }

    iterator end() && noexcept
    {
        assert(_journal != nullptr);
        return iterator{*_journal, size()};
    }
};

template <typename journal_t>
class journal_sequence<journal_t>::iterator
{
private:

    using journal_iterator = std::ranges::iterator_t<journal_t const &>;
    using journal_entry_t = std::iter_reference_t<journal_iterator>;
    using ref_sequence_t = std::tuple_element_t<1, std::remove_reference_t<journal_entry_t>>;
    using sequence_iterator = std::ranges::iterator_t<ref_sequence_t>;

    journal_iterator _dict_first{};
    journal_iterator _dict_last{};
    journal_iterator _dict_it{};
    size_t _position{}; // the current global position
    size_t _previous_switch{}; // the next position to switch in increment
    size_t _next_switch{}; // the next position to switch in increment
    sequence_iterator _sequence_it{}; // the current segment iterator

public:

    using value_type = std::iter_value_t<sequence_iterator>;
    using reference = std::iter_reference_t<sequence_iterator>;
    using pointer = std::add_pointer_t<value_type const>;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;

    iterator() = default;
    iterator(iterator const &) = default;
    iterator(iterator &&) = default;
    iterator & operator=(iterator const &) = default;
    iterator & operator=(iterator &&) = default;
    ~iterator() = default;

    iterator(journal_t const & journal, size_t const position) :
        _dict_first{journal.begin()},
        _dict_last{journal.end()},
        _dict_it{_dict_last},
        _position{position},
        _previous_switch{position},
        _next_switch{position}
    {
        // _entry_last = _host->_dictionary.end();
        if (position == 0) {
            _dict_it = _dict_first;
            _next_switch = journal_t::entry_last(*_dict_it);
            _sequence_it = std::ranges::begin(journal_t::entry_value(*_dict_it));
        }
    }

    iterator base() const & noexcept {
        return *this;
    }

    iterator && base() && noexcept {
        return std::move(*this);
    }

    /*!\name Element access
     * \{
     */
    reference operator*() const noexcept
    {
        return *_sequence_it;
    }

    pointer operator->() const noexcept
    {
        return std::addressof(this->operator*());
    }

    reference operator[](difference_type const offset) const noexcept
    {
        return *(*this + offset);
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    iterator & operator++() noexcept
    {
        ++_sequence_it; // now we may at the end of the segment
        if (++_position == _next_switch) [[unlikely]] { // we may or may not update the index
            if (++_dict_it != _dict_last) [[likely]]
                init_segment_begin();
        }

        return *this;
    }

    iterator operator++(int) noexcept(std::is_nothrow_copy_constructible_v<iterator>)
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    iterator & operator+=(difference_type offset) noexcept
    {
        // if offset is negative, use the preceding elements of current segment to check if binary search is needed.
        // if offset is positive use remaining elements of current segment to check if binary search is needed.
        constexpr auto element_proj = [] (auto && element) -> uint32_t { return journal_t::entry_last(element); };
        _position += offset;
        if (_position < _previous_switch || _next_switch <= _position) [[likely]] {
            if (_dict_it = std::ranges::upper_bound(_dict_first, _dict_last, _position, std::less{}, element_proj);
                _dict_it != _dict_last) {
                init_segment_begin();
                std::ranges::advance(_sequence_it, _position - journal_t::entry_first(*_dict_it));
            } else {
                _previous_switch = _position;
            }
        } else {
            _sequence_it += offset;
        }

        return *this;
    }

    iterator operator+(difference_type const offset) const noexcept
    {
        iterator tmp{*this};
        return tmp += offset;
    }

    friend iterator operator+(difference_type const offset, iterator const & rhs) noexcept
    {
        return rhs + offset;
    }

    iterator & operator--() noexcept
    {
        if (_position == _previous_switch) [[unlikely]] {
            --_dict_it;
            init_segment_end();
        }
        --_position;
        --_sequence_it;

        return *this;
    }

    iterator operator--(int) noexcept(std::is_nothrow_copy_constructible_v<iterator>)
    {
        iterator tmp{*this};
        --(*this);
        return tmp;
    }

    iterator & operator-=(difference_type offset) noexcept
    {
        return (*this) += -offset;
    }

    iterator operator-(difference_type offset) const noexcept
    {
        iterator tmp{*this};
        return tmp -= offset;
    }

    difference_type operator-(iterator const & rhs) const noexcept
    {
        return _position - rhs._position;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    bool operator==(iterator const & rhs) const noexcept
    {
        return _position == rhs._position;
    }

    auto operator<=>(iterator const & rhs) const noexcept
    {
        return _position <=> rhs._position;
    }
    //!\}

private:

    void init_segment_begin() noexcept {
        _previous_switch = journal_t::entry_first(*_dict_it);
        _next_switch = journal_t::entry_last(*_dict_it);
        _sequence_it = std::ranges::begin(journal_t::entry_value(*_dict_it));
    }

    void init_segment_end() noexcept {
        _previous_switch = journal_t::entry_first(*_dict_it);
        _next_switch = journal_t::entry_last(*_dict_it);
        _sequence_it = std::ranges::end(journal_t::entry_value(*_dict_it));
    }
};

}  // namespace libjst

namespace std::ranges {

template <typename journal_t>
inline constexpr bool enable_borrowed_range<libjst::journal_sequence<journal_t>> = true;

} // namespace std::ranges


