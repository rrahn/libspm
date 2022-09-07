// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an implementation for considering the jst as a range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>

#include <libjst/traversal/lazy_tree.hpp>

namespace libjst
{

    template <typename node_t>
    class jst_range
    {
    private:
        using branch_node_t = typename node_t::branch_node_type; // for the moment!
        using variant_branch_t = lazy_tree<branch_node_t>;

        class iterator;
        using sentinel = std::default_sentinel_t;

        variant_branch_t _current_branch{};
        node_t _current_base_node{};

    public:

        template <typename _node_t, observable_stack ...subscriber>
            requires (std::same_as<std::remove_cvref_t<_node_t>, node_t>)
        explicit jst_range(_node_t && root, subscriber_t & ...subscriber) :
            _current_branch{subscriber...},
            _current_base_node{(_node_t &&) root}
        {
        }

        iterator begin() noexcept
        {
            return iterator{this};
        }

        sentinel end() const noexcept
        {
            return std::default_sentinel;
        }
    };

    template <typename node_t, typename sink_t>
    class jst_range<node_t, sink_t>::iterator
    {
    private:
        using branch_iterator = std::ranges::iterator_t<variant_branch_t>;
        using node_value_t = typename node_t::value_type; // node carries some value type.

        friend sentinel;
        jst_range * _host{};
        node_value_t const *_current_value{};
        branch_iterator _branch_iter{};
    public:

        using value_type = node_value_t;
        using reference = value_type const &;
        using pointer = value_type const *;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;
        explicit iterator(jst_range * host) noexcept : _host{host}
        {
            if (_host != nullptr) {
                set_current_value(_host->_current_base_node);
            }
        }
        iterator(iterator const &) = delete;
        iterator(iterator &&) = default;
        iterator &operator=(iterator const &) = delete;
        iterator &operator=(iterator &&) = default;

        reference operator*() const noexcept
        {
            assert(_host != nullptr);
            return *_current_value;
        }

        pointer operator->() const noexcept
        {
            assert(_host != nullptr);
            return _current_value;
        }

        iterator & operator++() noexcept
        {
            assert(_host != nullptr);

            if (_branch_iter == std::ranges::end(_host->_current_branch)) {
                next_node();
            } else {
                if (++_branch_iter == std::ranges::end(_host->_current_branch)) {
                    set_current_value(_host->_current_base_node);
                else
                    set_current_value(*_branch_iter);
            }
            return *this;
        }

        node_value_t operator++(int) noexcept(std::is_nothrow_copy_constructible_v<node_value_t>)
        {
            assert(_host != nullptr);
            node_value_t tmp{*(*this)};
            ++(*this);
            return tmp;
        }

        constexpr bool operator==(std::default_sentinel_t const &) const noexcept
        {
            return _host != nullptr && !_host->_current_base_node.has_value();
        }

    private:
        template <typename current_node_t>
        void set_current_value(current_node_t const & node) noexcept
        {
            _current_value = std::addressof(*node);
        }

        void next_node() noexcept(std::is_nothrow_movable_v<branch_node_t>)
        {
            if (auto opt_branch = _host->_current_base_node.next(); opt_branch.has_value()) {
                _host->_current_branch.reset(std::move(*opt_branch)); // now we start a new variant branch cycle.
                _branch_iter = std::ranges::begin(_host->_current_branch);
                set_current_value(*_branch_iter);
            } else {
                set_current_value(_host->_current_base_node);
            }
        }
    };
}  // namespace libjst
