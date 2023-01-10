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

#include <concepts>
#include <ranges>

#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/take.hpp>

#include <libjst/journal_sequence.hpp>

namespace libjst::detail
{

template <std::ranges::viewable_range sequence_t>
using subrange_t = decltype(std::declval<sequence_t>() | seqan3::views::slice(0, 1));

} // namespace libjst::detail

namespace libjst {
// insertion sequence is different then reference sequence
template <std::integral position_t, std::ranges::forward_range sequence_t>
//!\cond
    requires std::ranges::viewable_range<sequence_t>
//!\endcond
class journal
{
private:

    using ref_sequence_type = detail::subrange_t<sequence_t>;
    using journal_entry_type = std::pair<position_t, ref_sequence_type>;
    using dictionary_type = std::vector<journal_entry_type>;
    using dictionary_iterator = std::ranges::iterator_t<dictionary_type>;

    template <bool is_const>
    class _iterator;

    template <typename>
    friend class journal_sequence;

    size_t _sequence_size{};
    dictionary_type _dictionary{};

public:

    using key_type = position_t;
    using mapped_type = ref_sequence_type;
    using value_type = std::pair<key_type const, detail::subrange_t<sequence_t const>>;
    using reference = value_type;
    using iterator = _iterator<false>;
    using const_iterator = _iterator<true>;
    using journaled_sequence_type = journal_sequence<journal>;

    journal() = default;
    explicit journal(sequence_t sequence) : _sequence_size{std::ranges::size(sequence)} {
        _dictionary.emplace_back(position_t{}, to_mapped(std::forward<sequence_t>(sequence)));

    }

    template <typename insert_sequence_t>
        requires std::constructible_from<mapped_type, detail::subrange_t<insert_sequence_t>>
    iterator record_insertion(key_type const position, insert_sequence_t && sequence)
    {
        if (std::ranges::empty(sequence))
            return iterator{find_entry(position)};

        assert(position <= static_cast<key_type>(sequence_size()));

        dictionary_iterator insert_it{};
        if (_dictionary.empty() || position == key_type{}) { // direct insertion
            auto sequence_size = std::ranges::size(sequence);
            insert_it = _dictionary.emplace(std::ranges::begin(_dictionary),
                                            position,
                                            to_mapped(std::forward<insert_sequence_t>(sequence)));
            rebalance_dictionary(std::ranges::next(insert_it), sequence_size);
        } else { // split entry
            insert_it = record_insertion_impl(find_entry(position),
                                              position,
                                              to_mapped(std::forward<insert_sequence_t>(sequence)));
        }
        return iterator{insert_it};
    }

    iterator record_deletion(key_type const position, size_t const count)
    {
        if (count == 0)
            return iterator{find_entry(position)};

        assert(check_valid_range(position, position + count));

        return iterator{record_deletion_impl(find_entry(position), position, count)};
    }

    template <typename insert_sequence_t>
        requires std::constructible_from<mapped_type, detail::subrange_t<insert_sequence_t>>
    iterator record_substitution(key_type const position, insert_sequence_t && sequence)
    {
        if (std::ranges::empty(sequence))
            return iterator{find_entry(position)};

        assert(check_valid_range(position, position + std::ranges::size(sequence)));

        return iterator{record_substitution_impl(find_entry(position),
                                                 position,
                                                 to_mapped(std::forward<insert_sequence_t>(sequence)))};
    }

    size_t size() const noexcept
    {
        return _dictionary.size();
    }

    bool empty() const noexcept
    {
        return _dictionary.empty();
    }

    iterator begin() noexcept
    {
        return iterator{std::ranges::begin(_dictionary)};
    }

    const_iterator begin() const noexcept
    {
        return const_iterator{std::ranges::begin(_dictionary)};
    }

    iterator end() noexcept
    {
        return iterator{std::ranges::end(_dictionary)};
    }

    const_iterator end() const noexcept
    {
        return const_iterator{std::ranges::end(_dictionary)};
    }

    constexpr journaled_sequence_type sequence() const noexcept {
        return journaled_sequence_type{*this};
    }

protected:

    constexpr size_t & sequence_size() noexcept {
        return _sequence_size;
    }

    constexpr size_t sequence_size() const noexcept {
        return _sequence_size;
    }

    constexpr dictionary_type & dictionary() noexcept {
        return _dictionary;
    }

    constexpr dictionary_type const & dictionary() const noexcept {
        return _dictionary;
    }

    dictionary_iterator record_insertion_impl(dictionary_iterator original_entry,
                                              key_type const position,
                                              mapped_type segment)
    {
        std::ptrdiff_t const insert_size = std::ranges::ssize(segment);
        dictionary_iterator dict_it = emplace_entry_hint(original_entry, position, std::move(segment));
        rebalance_dictionary(dict_it + 1, insert_size);
        return dict_it;
    }

    dictionary_iterator record_deletion_impl(dictionary_iterator original_entry,
                                             key_type const first,
                                             size_t const count)
    {
        dictionary_iterator dict_it = erase_range(original_entry, first, first + count);
        rebalance_dictionary(dict_it, -static_cast<std::ptrdiff_t>(count));
        return dict_it;
    }

    dictionary_iterator record_substitution_impl(dictionary_iterator original_entry,
                                                 key_type const position,
                                                 mapped_type segment)
    {
        size_t const segment_size = std::ranges::size(segment);
        key_type const last = position + segment_size;
        auto dict_it = _dictionary.emplace(erase_range(original_entry, position, last), position, std::move(segment));
        assert(check_consistent_segments());
        return dict_it;
    }

    constexpr dictionary_iterator find_entry(key_type const position) noexcept {
        assert(!_dictionary.empty());
        reserve_buffer();
        dictionary_iterator last = std::ranges::prev(std::ranges::end(_dictionary));
        if (position <= entry_first(*last)) {
            last = lower_bound(_dictionary, position);
        }
        return last;
    }

    void rebalance_dictionary(dictionary_iterator first, std::ptrdiff_t const offset) noexcept {
        std::ranges::for_each(first, std::ranges::end(_dictionary), [offset] (auto & entry) {
            entry_first(entry) += offset;
        });
        sequence_size() += offset;
        assert(check_consistent_segments());
    }

private:

    void reserve_buffer() {
        size_t current_capacity = _dictionary.capacity();
        _dictionary.reserve(current_capacity << (current_capacity == _dictionary.size()));
    }

    bool check_valid_range(key_type const first, key_type const last) const noexcept
    {
        return first < last && last <= static_cast<key_type>(sequence_size());
    }

    constexpr dictionary_iterator emplace_entry_hint(dictionary_iterator hint,
                                                     key_type const insert_position,
                                                     mapped_type && segment) noexcept {
        assert(hint != std::ranges::end(_dictionary));
        auto & affected_entry = *hint;

        assert(entry_first(affected_entry) < insert_position);
        assert(insert_position <= entry_last(affected_entry));

        // Create local entries for insertion and right entry covering suffix of split entry.
        key_type const split_position = insert_position - entry_first(affected_entry);
        journal_entry_type split_entry_right{insert_position,
                                             entry_value(affected_entry) | seqan3::views::drop(split_position)};
        journal_entry_type insert_entry{insert_position, std::move(segment)};

        // Update dictionary
        entry_value(affected_entry) = entry_value(affected_entry) | seqan3::views::take(split_position);
        if (!entry_value(split_entry_right).empty()) [[likely]] {
            return _dictionary.insert(++hint, {std::move(insert_entry), std::move(split_entry_right)});
        } else {
            return _dictionary.insert(++hint, std::move(insert_entry));
        }
    }

    constexpr dictionary_iterator erase_range(dictionary_iterator entry_left_it,
                                              key_type first,
                                              key_type last) noexcept {
        journal_entry_type & entry_left = *entry_left_it;
        key_type const prefix_last = first - entry_first(entry_left);

        // A: erase infix - requires split
        if ((prefix_last > 0) && (last < entry_last(entry_left))) [[likely]] {
            key_type const suffix_first = last - entry_first(entry_left);
            assert(prefix_last < suffix_first);
            assert(suffix_first < entry_last(entry_left));

            mapped_type suffix = entry_value(entry_left) | seqan3::views::drop(suffix_first);
            entry_value(entry_left) = entry_value(entry_left) | seqan3::views::take(prefix_last);
            assert(!entry_value(entry_left).empty());
            assert(!suffix.empty());
            return _dictionary.emplace(++entry_left_it, last, std::move(suffix));
        } else { // B: erease suffix of left entry and infix of right entry
            // 1. find right entry
            dictionary_iterator entry_right_it =
                lower_bound(std::ranges::subrange{entry_left_it, std::ranges::end(_dictionary)}, last);

            // 2. remove suffix of entry_left and cache prefix of entry right.
            bool keep_prefix_left = prefix_last > 0;
            bool erase_right = last == entry_last(*entry_right_it);

            key_type suffix_first = last - entry_first(*entry_right_it);
            mapped_type suffix_right = entry_value(*entry_right_it) | seqan3::views::drop(suffix_first);
            entry_value(entry_left) = entry_value(entry_left) | seqan3::views::take(prefix_last);

            // 3. erase inner entries
            entry_right_it = _dictionary.erase(entry_left_it + keep_prefix_left, entry_right_it + erase_right);
            // 4. store cached suffix to right entry.
            if (!erase_right) {
                assert(entry_right_it != std::ranges::end(_dictionary));
                entry_first(*entry_right_it) += suffix_first;
                entry_value(*entry_right_it) = std::move(suffix_right);
            }
            return entry_right_it;
        }
    }

    bool check_consistent_segments() const noexcept
    {
        bool is_consistent{true};
        key_type last_end{};
        std::ranges::for_each(_dictionary, [&] (auto const & entry)
        {
            if (entry_first(entry) != last_end)
                is_consistent = false;

            last_end += std::ranges::size(entry_value(entry));
        });
        return is_consistent;
    }

    template <typename insert_sequence_t>
    constexpr mapped_type to_mapped(insert_sequence_t && sequence) noexcept {
        return std::forward<insert_sequence_t>(sequence) | seqan3::views::take(std::ranges::size(sequence));
    }

    template <typename entry_t>
    static constexpr decltype(auto) entry_value(entry_t && entry) noexcept {
        return get<1>(std::forward<entry_t>(entry));
    }

    template <typename entry_t>
    static constexpr decltype(auto) entry_first(entry_t && entry) noexcept {
        return get<0>(std::forward<entry_t>(entry));
    }

    template <typename entry_t>
    static constexpr key_type entry_last(entry_t && entry) noexcept {
        return entry_first(entry) + std::ranges::size(entry_value(entry));
    }

    template <typename dictionary_t>
    static constexpr auto lower_bound(dictionary_t && dictionary, key_type const key) noexcept
        -> std::ranges::iterator_t<dictionary_t> {
        constexpr auto element_proj = [] (auto && element) -> key_type {
            return entry_last(std::forward<decltype(element)>(element));
        };
        return std::ranges::lower_bound(std::forward<dictionary_t>(dictionary), key, std::less{}, element_proj);
    }
};

template <typename sequence_t>
journal(sequence_t &&) -> journal<uint32_t, sequence_t>;

template <std::integral position_t, std::ranges::forward_range sequence_t>
//!\cond
    requires std::ranges::viewable_range<sequence_t>
//!\endcond
template <bool is_const>
class journal<position_t, sequence_t>::_iterator
{
private:
    using maybe_const_dictionary_type = std::conditional_t<is_const, dictionary_type const, dictionary_type>;
    using dictionary_iterator = std::ranges::iterator_t<maybe_const_dictionary_type>;

    template <bool>
    friend class _iterator;

    dictionary_iterator _iter;
public:

    using value_type = typename journal::value_type;
    using reference = typename journal::reference;
    using difference_type = std::ranges::range_difference_t<maybe_const_dictionary_type>;
    using pointer = std::conditional_t<is_const, typename dictionary_type::const_pointer, typename dictionary_type::pointer>;
    using iterator_category = std::random_access_iterator_tag;

    constexpr _iterator() = default;
    constexpr explicit _iterator(dictionary_iterator iter) noexcept : _iter{std::move(iter)}
    {}

    constexpr _iterator(_iterator<!is_const> iter) noexcept
        requires is_const
        : _iter{std::move(iter._iter)}
    {}

    dictionary_iterator base() const &
    {
        return _iter;
    }

    dictionary_iterator base() &&
    {
        return std::move(_iter);
    }

    constexpr reference operator*() const noexcept
    {
        return *_iter;
    }

    constexpr pointer operator->() const noexcept
    {
        return std::addressof(*_iter);
    }

    constexpr _iterator & operator++() noexcept
    {
        ++_iter;
        return *this;
    }

    constexpr _iterator operator++(int) noexcept
    {
        _iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr _iterator & operator+=(difference_type const offset) noexcept
    {
        _iter += offset;
        return *this;
    }

    constexpr _iterator operator+(difference_type const offset) const noexcept
    {
        _iterator tmp{*this};
        return tmp += offset;
    }

    friend constexpr _iterator operator+(difference_type const offset, _iterator const & rhs) noexcept
    {
        return rhs + offset;
    }

    constexpr _iterator & operator--() noexcept
    {
        --_iter;
        return *this;
    }

    constexpr _iterator operator--(int) noexcept
    {
        _iterator tmp{*this};
        --(*this);
        return tmp;
    }

    constexpr _iterator & operator-=(difference_type const offset) noexcept
    {
        _iter -= offset;
        return *this;
    }

    constexpr _iterator operator-(difference_type const offset) const noexcept
    {
        _iterator tmp{*this};
        return tmp -= offset;
    }

    template <bool other_const>
    constexpr difference_type operator-(_iterator<other_const> const & rhs) const noexcept
    {
        return _iter - rhs._iter;
    }

    template <bool other_const>
    constexpr bool operator==(_iterator<other_const> const & rhs) const noexcept
    {
        return _iter == rhs._iter;
    }

    template <bool other_const>
    constexpr auto operator<=>(_iterator<other_const> const & rhs) const noexcept
    {
        return _iter <=> rhs._iter;
    }
};

}  // namespace libjst
