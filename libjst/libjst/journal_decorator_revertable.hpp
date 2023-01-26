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
#include <set>
#include <span>
#include <stack>
#include <type_traits>

#include <seqan3/utility/views/slice.hpp>

#include <libjst/journal_entry.hpp>
#include <libjst/utility/sorted_vector.hpp>

namespace libjst
{

template <typename value_t>
class stack : public std::stack<value_t, std::vector<value_t>> {

    using base_t = std::stack<value_t, std::vector<value_t>>;
public:

    void reserve(size_t const new_capacity) {
        base_t::c.reserve(new_capacity);
    }
};

template <typename journal_decorator_t>
class journal_decorator_revertable : public journal_decorator_t
{
private:

    enum operation_kind {
        substitution,
        insertion,
        deletion
    };

    using typename journal_decorator_t::dictionary_type;
    using typename journal_decorator_t::dictionary_iterator;

    using history_element_t = std::tuple<size_t, typename journal_decorator_t::segment_type, operation_kind>;

    using stack_t = stack<history_element_t>;

    stack_t _history{};

public:

    using typename journal_decorator_t::value_type;
    using typename journal_decorator_t::size_type;
    using typename journal_decorator_t::reference;
    using typename journal_decorator_t::const_reference;
    using typename journal_decorator_t::segment_type;

    journal_decorator_revertable() = delete;

    explicit journal_decorator_revertable(journal_decorator_t && base) : journal_decorator_t{std::move(base)}
    {
        _history.reserve(100);
    }

    constexpr reference operator[](std::ptrdiff_t const pos) noexcept {
        return begin()[pos];
    }

    constexpr const_reference operator[](std::ptrdiff_t const pos) const noexcept {
        return begin()[pos];
    }

    bool record_insertion(size_type const position, segment_type segment) {
        assert(!dictionary().empty());
        assert(position <= journal_decorator_t::size());
        // but we need the span of the original record to do the stuff.
        auto dict_it = journal_decorator_t::find_entry(position);
        _history.emplace(dict_it - journal_decorator_t::dictionary().begin(), (*dict_it).segment(), operation_kind::insertion);
        journal_decorator_t::record_insertion_impl(dict_it, position, std::move(segment));
        return true;
    }

    // TODO: first plus count!
    bool record_deletion(size_type const first_position, size_type const last_position) {
        assert(journal_decorator_t::check_valid_range(first_position, last_position));

        auto dict_it = journal_decorator_t::find_entry(first_position);
        _history.emplace(dict_it - journal_decorator_t::dictionary().begin(), (*dict_it).segment(), operation_kind::deletion);
        journal_decorator_t::record_deletion_impl(dict_it, first_position, last_position);
        return true;
    }

    bool record_substitution(size_type const position, segment_type segment)
    {
        assert(journal_decorator_t::check_valid_range(position, position + std::ranges::size(segment)));

        auto dict_it = journal_decorator_t::find_entry(position);
        _history.emplace(dict_it - journal_decorator_t::dictionary().begin(), (*dict_it).segment(), operation_kind::substitution);
        journal_decorator_t::record_substitution_impl(dict_it, position, std::move(segment));
        return true;
    }

    typename journal_decorator_t::iterator<true> begin() noexcept
    {
        return journal_decorator_t::begin();
    }

    typename journal_decorator_t::iterator<false> begin() const noexcept
    {
        return journal_decorator_t::begin();
    }

    typename journal_decorator_t::iterator<true> end() noexcept
    {
        return journal_decorator_t::end();
    }

    typename journal_decorator_t::iterator<false> end() const noexcept
    {
        return journal_decorator_t::end();
    }

    // reverts the last recorded operation.
    void revert() {
        assert(!_history.empty());
        revert_impl();
    }

private:

    void revert_impl() {
        auto [original_index, original_segment, operation] = std::move(_history.top());
        _history.pop();

        auto dict_it = journal_decorator_t::dictionary()._elements.begin() + original_index;
        (*dict_it).segment() = std::move(original_segment);

        ++dict_it;
        std::ptrdiff_t effective_size{};
        if (operation == operation_kind::deletion) {
            dict_it = dictionary_iterator{journal_decorator_t::dictionary()._elements.erase(dict_it)}; // remove next
            effective_size = std::ranges::ssize(original_segment);
        } else {
            dict_it = dictionary_iterator{journal_decorator_t::dictionary()._elements.erase(
                    dict_it, std::ranges::next(dict_it, 2, journal_decorator_t::dictionary()._elements.end()))};
            effective_size = (operation == operation_kind::substitution) * -std::ranges::ssize(original_segment);
        }

        journal_decorator_t::rebalance_dictionary(dict_it, effective_size);
    }

};

}  // namespace libjst
