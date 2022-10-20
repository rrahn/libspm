// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides path.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>

#include <libjst/structure/concept_jst.hpp>
#include <libjst/variant/variant_proxy_offset.hpp>
#include <libjst/journal.hpp>

namespace libjst
{

    template <typename source_sequence_t, typename variant_store_t>
    class journaled_path
    {
    private:
        using variant_iterator_t = std::ranges::iterator_t<variant_store_t>;
        using variant_t = std::iter_value_t<variant_iterator_t>;
        using position_t = variant_position_t<variant_t>;
        using journal_t = journal<position_t, source_sequence_t>;
        using offset_t = std::ptrdiff_t;

        class iterator;

        journal_t _journal{};
        variant_iterator_t _root{};
        variant_iterator_t _sink{};
        offset_t _offset{};

        journaled_path(journal_t journal, variant_iterator_t root, variant_iterator_t sink, offset_t offset) :
            _journal{std::move(journal)},
            _root{std::move(root)},
            _sink{std::move(sink)},
            _offset{std::move(offset)}
        {
            // find the next root iterator.
        }

    public:

        template <typename _source_sequence_t, typename _variant_store_t>
            requires (std::constructible_from<journal_t, _source_sequence_t> &&
                      std::constructible_from<variant_iterator_t, std::ranges::iterator_t<_variant_store_t>>)
        journaled_path(_source_sequence_t && source, _variant_store_t && store) :
            _journal{source},
            _root{std::ranges::begin(store)},
            _sink{std::ranges::end(store)}
        {}

        auto sequence() const noexcept
        {
            return _journal.sequence();
        }

        journaled_path alternative_path(iterator it) const
        {
            auto const & variant = *it;
            // copy and add variant to journal.
            journal_t alt_journal{_journal};
            if (libjst::is_insertion(variant)) {
                alt_journal.record_insertion(libjst::position(variant), libjst::insertion(variant));
            } else if (libjst::is_deletion(variant)) {
                alt_journal.record_deletion(libjst::position(variant), libjst::deletion(variant));
            } else {
                assert(libjst::is_replacement(variant));
                alt_journal.record_substitution(libjst::position(variant), libjst::insertion(variant));
            }
            offset_t effective_size = std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
            return journaled_path{std::move(alt_journal),
                                  find_next(std::move(it).base()),
                                  _sink,
                                  _offset + effective_size};
        }

        iterator begin() const noexcept
        {
            return iterator{_root, _offset};
        }

        iterator end() const noexcept
        {
            return iterator{_sink, _offset};
        }

    private:

        variant_iterator_t find_next(variant_iterator_t it) const
        {
            position_t const branch_position = libjst::position(*it);
            position_t const branch_end = branch_position + libjst::deletion(*it);
            it = std::ranges::find_if(++it, _sink, [&] (auto &&variant) {
                return !libjst::is_insertion(variant) || libjst::position(variant) != branch_position;
            });

            // second: if next variant is not already the next valid we need to search it.
            if (it != _sink && branch_end > libjst::position(*it)) {
                it = std::ranges::lower_bound(it, _sink, branch_end, std::less<>{},
                                                         [] (auto &&variant) { return libjst::position(variant); });
            }

            return it;
        }
    };

    template <typename source_sequence_t, typename variant_store_t>
    journaled_path(source_sequence_t &&, variant_store_t &&)
        -> journaled_path<std::remove_reference_t<source_sequence_t>, std::remove_reference_t<variant_store_t>>;

    template <typename source_sequence_t, typename variant_iterator_t>
    class journaled_path<source_sequence_t, variant_iterator_t>::iterator
    {
    private:

        variant_iterator_t _it{};
        offset_t _offset{};
    public:

        using value_type = offset_variant<std::iter_value_t<variant_iterator_t>>;
        using reference = value_type;
        using difference_type = std::iter_difference_t<variant_iterator_t>;
        using pointer = void;
        using iterator_category = std::random_access_iterator_tag; // maybe that of the variant iterator?

        constexpr iterator() = default;
        constexpr explicit iterator(variant_iterator_t it, offset_t offset) :
            _it{std::move(it)},
            _offset{std::move(offset)}
        {
        }

        constexpr variant_iterator_t base() const & noexcept
        {
            return _it;
        }

        constexpr variant_iterator_t base() && noexcept
        {
            return std::move(_it);
        }

        constexpr reference operator*() const noexcept
        {
            return offset_variant{*_it, _offset};
        }

        constexpr iterator & operator++() noexcept
        {
            ++_it;
            return *this;
        }

        constexpr iterator operator++(int) noexcept
        {
            iterator tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator & operator+=(difference_type offset) noexcept
        {
            _it += offset;
            return *this;
        }

        constexpr iterator operator+(difference_type offset) const noexcept
        {
            iterator tmp{*this};
            return tmp += offset;
        }

        constexpr friend iterator operator+(difference_type offset, iterator rhs) noexcept
        {
            return rhs += offset;
        }

        constexpr iterator & operator--() noexcept
        {
            --_it;
            return *this;
        }

        constexpr iterator operator--(int) noexcept
        {
            iterator tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator & operator-=(difference_type offset) noexcept
        {
            _it -= offset;
            return *this;
        }

        constexpr iterator operator-(difference_type offset) const noexcept
        {
            iterator tmp{*this};
            return tmp -= offset;
        }

        constexpr difference_type operator-(iterator const & rhs) noexcept
        {
            return _it - rhs._it;
        }

        constexpr bool operator==(iterator const & rhs) const noexcept
        {
            return _it == rhs._it;
        }

        constexpr std::strong_ordering operator<=>(iterator const & rhs) const noexcept
        {
            return _it <=> rhs._it;
        }
    };
}  // namespace libjst
