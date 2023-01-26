// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>
#include <set>
#include <span>
#include <type_traits>

#include <seqan3/utility/views/slice.hpp>

#include <libjst/journal_entry.hpp>
#include <libjst/utility/sorted_vector.hpp>

namespace libjst
{

// TODO: rename to journaled_sequence with sequence interface.
 // uses internally the journal, which can be specified with some other journal type if we want to.
template <std::ranges::view segment_t>
class journal_decorator
{
protected:

    using journal_entry_type = detail::journal_entry<segment_t>;
    using dictionary_type = sorted_vector<journal_entry_type, std::less<void>>;
    using dictionary_iterator = typename dictionary_type::iterator;
    using key_compare = typename dictionary_type::key_compare;

private:
    dictionary_type _dictionary{};
    size_t _size{};

public:

    template <bool is_const_range>
    class iterator;

    using value_type = std::ranges::range_value_t<segment_t>;
    using size_type = typename journal_entry_type::size_type;
    using reference = std::ranges::range_reference_t<segment_t>;
    using const_reference = std::ranges::range_reference_t<segment_t const>;
    using segment_type = segment_t;

    journal_decorator() = default;
    journal_decorator(journal_decorator const &) = default;
    journal_decorator(journal_decorator &&) = default;
    journal_decorator & operator=(journal_decorator const &) = default;
    journal_decorator & operator=(journal_decorator &&) = default;
    ~journal_decorator() = default;

    explicit journal_decorator(segment_type initial_segment, size_t const initial_cacpacity = 32) noexcept
    {
        _dictionary._elements.reserve(initial_cacpacity);
        _size = std::ranges::size(initial_segment);
        if (_size > 0) // only add entry if initial segment is not empty.
            _dictionary.emplace(0, std::move(initial_segment));
    }

    constexpr reference operator[](std::ptrdiff_t const pos) noexcept {
        return begin()[pos];
    }

    constexpr const_reference operator[](std::ptrdiff_t const pos) const noexcept {
        return begin()[pos];
    }

    bool record_insertion(size_type const position, segment_type segment)
    {
        assert(position <= size());

        if (_dictionary.empty() || position == 0) { // direct insertion
            auto segment_size = std::ranges::size(segment);
            _dictionary._elements.emplace(_dictionary._elements.begin(), position, std::move(segment));
            rebalance_dictionary(std::ranges::next(_dictionary.begin()), segment_size);
        } else { // split entry
            record_insertion_impl(find_entry(position), position, std::move(segment));
        }
        return true;
    }

    bool record_deletion(size_type const first_position, size_type const last_position)
    {
        assert(check_valid_range(first_position, last_position));

        record_deletion_impl(find_entry(first_position), first_position, last_position);
        return true;
    }

    bool record_substitution(size_type const position, segment_type segment)
    {
        assert(check_valid_range(position, position + std::ranges::size(segment)));

        record_substitution_impl(find_entry(position), position, std::move(segment));
        return true;
    }

    // bool record_substitution_back(size_type const position, segment_type segment)
    // {
    //     reserve_buffer();
    //     size_type const segment_size = std::ranges::size(segment);
    //     size_type const last_position = position + segment_size;

    //     assert(check_valid_range(position, last_position));

    //     auto & last_entry = _dictionary._elements.back();

    //     // either inside of entry or at end or (at first in case of inserting at position 0).
    //     size_type const split_position = position - last_entry.segment_begin_position();
    //     size_type const segment_end_position = last_entry.segment_end_position();

    //     // Create local entries for insertion and right entry covering segment suffix of split entry.
    //     // journal_entry_type split_entry_right{last_position, last_entry.segment() |
    //     //                                         seqan3::views::slice(split_position + segment_size, segment_end_position)};
    //     // journal_entry_type insert_entry{split_position, std::move(segment)};

    //     // Update the entries in the dictionary!
    //     // we need to insert two at once:
    //     segment_type suffix = last_entry.segment()
    //                         | seqan3::views::slice(split_position + segment_size, segment_end_position);
    //     last_entry.segment() = last_entry.segment() | seqan3::views::slice(0, split_position);
    //     _dictionary._elements.emplace_back(position, std::move(segment));
    //     _dictionary._elements.emplace_back(last_position, std::move(suffix)); // | seqan3::views::slice(split_position + segment_size, segment_end_position));
    //      // ;
    //     // return dictionary_iterator{_dictionary._elements.insert(++hint.base(), {std::move(insert_entry), std::move(split_entry_right)})};
    //     // if (split_position != segment_end_position) {
    //     // } else {
    //     //     return dictionary_iterator{_dictionary._elements.insert(++hint.base(), std::move(insert_entry))};
    //     // }


    //     // _dictionary._elements.emplace(erase_range(position, last_position).base(), position, segment);

    //     assert(check_consistent_segments());
    //     return true;
    // }

    auto insert(iterator<true> const position,
                std::ranges::iterator_t<segment_t> first,
                std::ranges::iterator_t<segment_t> last) {
        return record_insertion(position - begin(), segment_type{first, last});
    }

    auto erase(iterator<true> const first, iterator<true> const last) {
        return record_deletion(first - begin(), last - begin());
    }

    // auto replace(size_type const position, [[maybe_unused]] size_t const count, segment_type segment) {
    //     assert(count == segment.size());
    //     return record_substitution(position, segment);
    // }

    size_type size() const noexcept
    {
        return _size;
    }

    bool empty() const noexcept
    {
        return _dictionary.empty();
    }

    iterator<true> begin() noexcept
    {
        return iterator<true>{_dictionary, 0};
    }

    iterator<true> begin() const noexcept
    {
        return iterator<true>{_dictionary, 0};
    }

    iterator<true> end() noexcept
    {
        return iterator<true>{_dictionary, size()};
    }

    iterator<true> end() const noexcept
    {
        return iterator<true>{_dictionary, size()};
    }

protected:

    constexpr dictionary_type & dictionary() noexcept {
        return _dictionary;
    }

    constexpr dictionary_type const & dictionary() const noexcept {
        return _dictionary;
    }

    dictionary_iterator record_insertion_impl(dictionary_iterator original_entry,
                                              size_type const position,
                                              segment_type segment)
    {
        std::ptrdiff_t const insert_size = std::ranges::ssize(segment);
        dictionary_iterator dict_it = emplace_entry_hint(original_entry, position, std::move(segment));
        rebalance_dictionary(dict_it + 1, insert_size);
        return dict_it;
    }

    dictionary_iterator record_deletion_impl(dictionary_iterator original_entry,
                                             size_type const first_position,
                                             size_type const last_position)
    {
        std::ptrdiff_t deletion_size = -static_cast<std::ptrdiff_t>(last_position - first_position);
        dictionary_iterator dict_it = erase_range(original_entry, first_position, last_position);
        rebalance_dictionary(dict_it, deletion_size);
        return dict_it;
    }

    dictionary_iterator record_substitution_impl(dictionary_iterator original_entry,
                                                 size_type const position,
                                                 segment_type segment)
    {
        size_type const segment_size = std::ranges::size(segment);
        size_type const last_position = position + segment_size;

        auto dict_it = _dictionary._elements.emplace(erase_range(original_entry, position, last_position).base(),
                                                     position,
                                                     std::move(segment));

        assert(check_consistent_segments());
        return dictionary_iterator{std::move(dict_it)};
    }

    constexpr dictionary_iterator find_entry(size_type const position) noexcept {
        assert(!_dictionary._elements.empty());
        reserve_buffer();
        dictionary_iterator last = std::ranges::prev(_dictionary.end());
        if (position <= (*last).segment_begin_position()) {
            last = _dictionary.lower_bound(position);
        }
        return last;
    }

    void rebalance_dictionary(dictionary_iterator first, std::ptrdiff_t const offset) noexcept {
        std::ranges::for_each(first.base(), _dictionary._elements.end(), [offset] (journal_entry_type & entry) {
            entry.segment_begin_position() += offset;
        });
        _size += offset;
        assert(check_consistent_segments());
    }
private:

    void reserve_buffer() {
        size_t current_capacity = _dictionary._elements.capacity();
        _dictionary._elements.reserve(current_capacity << (current_capacity == _dictionary._elements.size()));
    }

    bool check_valid_range(size_type const first_position, size_type const last_position) const noexcept
    {
        return first_position < last_position && last_position <= size();
    }

    void update_size(int32_t const event_size) noexcept
    {
        _size += event_size;
    }

    constexpr dictionary_iterator emplace_entry_hint(dictionary_iterator hint,
                                                     size_type const insert_position,
                                                     segment_type && segment) noexcept {
        assert(hint != _dictionary.end());
        auto & affected_entry = *hint;

        assert(affected_entry.segment_begin_position() < insert_position);
        assert(insert_position <= affected_entry.segment_end_position());
        // either inside of entry or at end or (at first in case of inserting at position 0).
        size_type const split_position = insert_position - affected_entry.segment_begin_position();
        // size_type const segment_end_position = affected_entry.segment_end_position();

        // Create local entries for insertion and right entry covering segment suffix of split entry.
        journal_entry_type split_entry_right{insert_position,
                                             affected_entry.segment() | std::views::drop(split_position)};
        journal_entry_type insert_entry{insert_position, std::move(segment)};

        // Update the entries in the dictionary!
        affected_entry.segment() = affected_entry.segment() | std::views::take(split_position);
        // we need to insert two at once:
        if (!split_entry_right.segment().empty()) [[likely]] {
            return dictionary_iterator{_dictionary._elements.insert(++hint.base(), {std::move(insert_entry), std::move(split_entry_right)})};
        } else {
            return dictionary_iterator{_dictionary._elements.insert(++hint.base(), std::move(insert_entry))};
        }
    }

    constexpr dictionary_iterator erase_range(dictionary_iterator entry_left_it,
                                              size_type first_position,
                                              size_type last_position) noexcept {
        journal_entry_type & entry_left = *entry_left_it;
        size_type const prefix_position = first_position - entry_left.segment_begin_position();

        // A: erase infix - requires split
        if ((prefix_position > 0) && (last_position < entry_left.segment_end_position())) [[likely]] {
            size_type const suffix_begin_position = last_position - entry_left.segment_begin_position();
            assert(prefix_position < suffix_begin_position);
            assert(suffix_begin_position < entry_left.segment_end_position());

            segment_type suffix = entry_left.segment() | seqan3::views::slice(suffix_begin_position,
                                                                                  entry_left.segment_end_position());
            entry_left.segment() = entry_left.segment() | seqan3::views::slice(0, prefix_position);
            assert(!entry_left.segment().empty());
            assert(!suffix.empty());
            return dictionary_iterator{_dictionary._elements.emplace(++entry_left_it.base(),
                                                                     last_position,
                                                                     std::move(suffix))};
        } else { // B: erease suffix of left entry and infix of right entry
            // 1. find right entry
            auto entry_right_it = std::ranges::lower_bound(entry_left_it.base(),
                                                           _dictionary._elements.end(),
                                                           last_position,
                                                           key_compare{});

            // 2. remove suffix of entry_left and cache prefix of entry right.
            bool keep_prefix_left = prefix_position > 0;
            bool erase_right = last_position == (*entry_right_it).segment_end_position();

            size_type suffix_position = last_position - (*entry_right_it).segment_begin_position();
            segment_type suffix_right = (*entry_right_it).segment()
                                      | seqan3::views::slice(suffix_position,
                                                             (*entry_right_it).segment_end_position());
            entry_left.segment() = entry_left.segment() | seqan3::views::slice(0, prefix_position);

            // 3. erase inner entries
            entry_right_it = _dictionary._elements.erase(entry_left_it.base() + keep_prefix_left,
                                                         entry_right_it + erase_right);
            // 4. store cached suffix to right entry.
            if (!erase_right) {
                assert(entry_right_it != _dictionary._elements.end());
                (*entry_right_it).segment_begin_position() += suffix_position;
                (*entry_right_it).segment() = std::move(suffix_right);
            }
            return dictionary_iterator{entry_right_it};
        }
    }

    bool check_consistent_segments() const noexcept
    {
        bool is_consistent{true};
        size_type last_end = 0;
        std::ranges::for_each(_dictionary, [&] (auto const & entry)
        {
            if (entry.segment_begin_position() != last_end)
                is_consistent = false;

            last_end += entry.segment_size();
        });
        return is_consistent;
    }
};

template <std::ranges::contiguous_range range_t>
journal_decorator(range_t &&) -> journal_decorator<std::remove_reference_t<range_t>>;

template <std::ranges::view segment_t>
template <bool is_const_range>
class journal_decorator<segment_t>::iterator
{
private:

    using maybe_const_segment_t = std::conditional_t<is_const_range, segment_t const, segment_t>;
    using maybe_const_dictionary_t = std::conditional_t<is_const_range, dictionary_type const, dictionary_type>;

    using segment_iterator = std::ranges::iterator_t<maybe_const_segment_t>;
    using dictionary_iterator = std::ranges::iterator_t<maybe_const_dictionary_t>;
    using journal_entry_t = std::iter_value_t<dictionary_iterator>;
    using segment_diff_t = std::iter_difference_t<segment_iterator>;

    maybe_const_dictionary_t * _dictionary{};

    dictionary_iterator _dict_it{};  // the iterator to the dictionary
    size_t _position{}; // the current global position
    size_t _previous_switch{}; // the next position to switch in increment
    size_t _next_switch{}; // the next position to switch in increment
    segment_iterator _sequence_it{}; // the current segment iterator

    template <bool>
    friend class iterator;

public:

    using value_type = std::ranges::range_value_t<maybe_const_segment_t>;
    using reference = std::ranges::range_reference_t<maybe_const_segment_t>;
    using pointer = std::add_pointer_t<value_type>;
    using difference_type = std::ranges::range_difference_t<maybe_const_segment_t>;
    using iterator_category = std::random_access_iterator_tag; // we should not but well.

    iterator() = default;
    iterator(iterator const &) = default;
    iterator(iterator &&) = default;
    iterator & operator=(iterator const &) = default;
    iterator & operator=(iterator &&) = default;
    ~iterator() = default;

    iterator(maybe_const_dictionary_t & dictionary, size_t const position) :
        _dictionary{std::addressof(dictionary)},
        _dict_it{std::ranges::end(*_dictionary)},
        _position{position},
        _previous_switch{position},
        _next_switch{position}
    {
        // _entry_last = _host->_dictionary.end();
        if (position == 0) {
            _dict_it = std::ranges::begin(*_dictionary);
            _next_switch = (*_dict_it).segment_end_position();
            _sequence_it = std::ranges::begin((*_dict_it).segment());
        }

        // if (!_host->_dictionary.empty()) {
        //     // _entry_last = std::ranges::prev(_host->_dictionary.end());
        //     if (is_end) {
        //         _entry_it = _host->_dictionary.end();
        //         _segment_position = _host->size();
        //     } else {
        //         _entry_it = _host->_dictionary.begin();
        //         _segment_position = 0;
        //     }
        // }
    }

    iterator(iterator<!is_const_range> other)
        requires is_const_range
    :
        _dictionary{other._dictionary},
        _dict_it{other._dict_it},
        _position{other._position},
        _previous_switch{other._previous_switch},
        _next_switch{other._next_switch},
        _sequence_it{other._sequence_it}
    {}

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
            if (++_dict_it != std::ranges::end(*_dictionary)) [[likely]]
                init_segment_begin();
        }
        // if (_segment_position == _next_switch) {
        //     ++_entry_it;
        // }
        // _entry_it = std::ranges::next(_entry_it, ++_segment_position == _entry_it->segment_end_position(), _entry_last);
        // _entry_it = (++_segment_position == _entry_it->segment_end_position()) ?
        // if (++_segment_position == _entry_it->segment_end_position()) {
        //     ++_entry_it;
        // }

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
        _position += offset;
        if (_position < _previous_switch || _next_switch <= _position) [[likely]] {
            if (_dict_it = _dictionary->upper_bound(_position); _dict_it != std::ranges::end(*_dictionary)) {
                init_segment_begin();
                std::ranges::advance(_sequence_it, _position - (*_dict_it).segment_begin_position());
            } else {
                _previous_switch = _position;
            }
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

    template <bool other_const_range>
    difference_type operator-(iterator<other_const_range> const & rhs) const noexcept
    {
        return _position - rhs._position;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    template <bool other_const_range>
    bool operator==(iterator<other_const_range> const & rhs) const noexcept
    {
        return _position == rhs._position;
    }

    template <bool other_const_range>
    auto operator<=>(iterator<other_const_range> const & rhs) const noexcept
    {
        return _position <=> rhs._position;
    }
    //!\}

private:

    void init_segment_begin() noexcept {
        _previous_switch = (*_dict_it).segment_begin_position();
        _next_switch = (*_dict_it).segment_end_position();
        _sequence_it = std::ranges::begin((*_dict_it).segment());
    }

    void init_segment_end() noexcept {
        _previous_switch = (*_dict_it).segment_begin_position();
        _next_switch = (*_dict_it).segment_end_position();
        _sequence_it = std::ranges::end((*_dict_it).segment());
    }
};

}  // namespace libjst
