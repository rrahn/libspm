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
#include <memory>
#include <ranges>

#include <seqan3/range/views/slice.hpp>

#include <libjst/container/concept_jst.hpp>
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

        class node_type;
        class iterator;
        using sentinel = std::default_sentinel_t;

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
        }

    public:

        template <typename _source_sequence_t, typename _variant_store_t>
            requires (std::constructible_from<journal_t, _source_sequence_t> &&
                      std::constructible_from<variant_iterator_t, std::ranges::iterator_t<_variant_store_t>>)
        journaled_path(_source_sequence_t const & source, _variant_store_t && store) :
            _journal{source},
            _root{std::ranges::begin(store)},
            _sink{std::ranges::end(store)}
        {}

        auto label() const noexcept
        {
            return _journal.sequence();
        }

        journaled_path alternative_path(node_type const &node) const
        {
            variant_t const & variant = *node._next_variant; // now we get a node!
            position_t const journal_position = libjst::position(variant) + _offset;
            // copy and add variant to journal.
            journal_t alt_journal{_journal};
            if (libjst::is_insertion(variant)) {
                alt_journal.record_insertion(journal_position, libjst::insertion(variant));
            } else if (libjst::is_deletion(variant)) {
                alt_journal.record_deletion(journal_position, libjst::deletion(variant));
            } else {
                assert(libjst::is_replacement(variant));
                alt_journal.record_substitution(journal_position, libjst::insertion(variant));
            }
            offset_t effective_size = std::ranges::size(libjst::insertion(variant)) - libjst::deletion(variant);
            return journaled_path{std::move(alt_journal),
                                  find_next(node._next_variant),
                                  _sink,
                                  _offset + effective_size};
        }

        iterator begin() const & noexcept
        {
            return iterator{node_type{*this, 0}};
        }

        iterator begin() && noexcept
        {
            return iterator{node_type{std::move(*this), 0}};
        }

        sentinel end() const & noexcept
        {
            return {};
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

    template <typename source_sequence_t, typename variant_store_t>
    class journaled_path<source_sequence_t, variant_store_t>::node_type
    {
    private:

        friend journaled_path;

        std::shared_ptr<journaled_path> _self_managed_path{};
        journaled_path const * _path{};

        variant_iterator_t _next_variant{};
        position_t _label_begin_position{}; // thus the end node has the same path and does not need to be allocated twice.

        constexpr node_type() = default;
        // construction to external object.
        explicit constexpr node_type(journaled_path const & path, position_t position) :
            _path{std::addressof(path)},
            _next_variant{_path->_root},
            _label_begin_position(position)
        {
        }

        explicit constexpr node_type(journaled_path && path, position_t position) :
            _self_managed_path{std::make_shared<journaled_path>(std::move(path))},
            _path{_self_managed_path.get()},
            _next_variant{_path->_root},
            _label_begin_position(position)
        {}

    public:
        constexpr auto label() const noexcept
        {
            return path() | seqan3::views::slice(begin_position(), label_end_position());
        }

        constexpr auto label_size() const noexcept
        {
            return label_end_position() - begin_position();
        }

        // start position inside the path
        // sink starts with position size(path)
        constexpr position_t begin_position() const noexcept
        {
            return _label_begin_position;
        }

        constexpr position_t end_position() const noexcept
        {
            return label_end_position();
        }

        constexpr node_type & next() noexcept
        {
            _label_begin_position = label_end_position();
            ++_next_variant;
            return *this;
        }

        constexpr node_type alt() const
        {
            assert(!is_leaf());
            return node_type{_path->alternative_path(*this), label_end_position()};
        }

        constexpr variant_iterator_t next_variant() const noexcept
        {
            return _next_variant;
        }

        constexpr bool is_leaf() const noexcept
        {
            return _next_variant == path()._sink && begin_position() != label_end_position();
        }

        constexpr operator bool() const noexcept
        {
            return _next_variant != path()._sink || begin_position() != label_end_position();
        }

        journaled_path const & path() const noexcept
        {
            assert(_path != nullptr);
            return *_path;
        }

    private:

        constexpr position_t label_end_position() const noexcept
        {
            return (_next_variant != path()._sink) ? libjst::position(*_next_variant) + path()._offset
                                                   : std::ranges::ssize(path().label());
        }
    };

    template <typename source_sequence_t, typename variant_iterator_t>
    class journaled_path<source_sequence_t, variant_iterator_t>::iterator
    {
    private:

        node_type _node{};
    public:

        using value_type = node_type;
        using reference = node_type const &;
        using difference_type = std::iter_difference_t<variant_iterator_t>;
        using pointer = node_type const *;
        using iterator_category = std::input_iterator_tag; // maybe that of the variant iterator?

        constexpr iterator() = default;
        constexpr explicit iterator(node_type node) : _node{std::move(node)}
        {
        }

        constexpr reference operator*() const noexcept
        {
            return _node;
        }

        constexpr iterator & operator++() noexcept
        {
            _node.next();
            return *this;
        }

        constexpr iterator operator++(int) noexcept
        {
            iterator tmp{*this};
            ++(*this);
            return tmp;
        }

        // constexpr iterator & operator+=(difference_type offset) noexcept
        // {
        //     _it += offset;
        //     return *this;
        // }

        // constexpr iterator operator+(difference_type offset) const noexcept
        // {
        //     iterator tmp{*this};
        //     return tmp += offset;
        // }

        // constexpr friend iterator operator+(difference_type offset, iterator rhs) noexcept
        // {
        //     return rhs += offset;
        // }

        // constexpr iterator & operator--() noexcept
        // {
        //     --_it;
        //     return *this;
        // }

        // constexpr iterator operator--(int) noexcept
        // {
        //     iterator tmp{*this};
        //     --(*this);
        //     return tmp;
        // }

        // constexpr iterator & operator-=(difference_type offset) noexcept
        // {
        //     _it -= offset;
        //     return *this;
        // }

        // constexpr iterator operator-(difference_type offset) const noexcept
        // {
        //     iterator tmp{*this};
        //     return tmp -= offset;
        // }

        // constexpr difference_type operator-(iterator const & rhs) noexcept
        // {
        //     return _it - rhs._it;
        // }

        constexpr bool operator==(sentinel const &) const noexcept
        {
            return !_node;
        }

        // constexpr std::strong_ordering operator<=>(iterator const & rhs) const noexcept
        // {
        //     return _it <=> rhs._it;
        // }
    };
}  // namespace libjst

namespace std::ranges {

template <typename source_sequence_t, typename variant_store_t>
inline constexpr bool enable_borrowed_range<libjst::journaled_path<source_sequence_t, variant_store_t>> = true;

} // namespace std::ranges



