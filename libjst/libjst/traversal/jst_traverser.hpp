// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides jst base traversal.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

    template <typename jst_t>
    class jst_traverser // implementing dfs_traversal?
    {
    private:

        using jst_node_type = typename jst_t::node_type;

        class context;
        class iterator;

        jst_t _jst{};
        std::stack<jst_node_type> _active_branch{};

    public:
        constexpr explicit jst_traverser(jst_t jst) noexcept : _jst{std::move(jst)}
        {
            _active_branch.push_back(_jst.sink()); // add the sink as last node.
            _active_branch.push_back(_jst.root()); // add the root as starting node.
        }

        iterator begin() noexcept {
            return iterator{*this};
        }

        std::default_sentinel_t end() noexcept {
            return iterator{*this};
        }
    };

    template <typename jst_t>
    class jst_traverser<jst_t>::context : public traversal_context {
    private:

    public:

        context() = default;

        // what can be added to this?
    };

    template <typename jst_t>
    class jst_traverser<jst_t>::iterator
    {
    private:

        jst_traverser * _host{};

        constexpr explicit iterator(jst_traverser & host) noexcept : _host{host} {
            ++(*this); // move to the first active node.
        }

    public:

        using value_type = context;
        using reference = context;
        using pointer = void;
        using diference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        iterator() = default;
        iterator(iterator const &) = delete;
        iterator(iterator &&) = default;
        iterator & operator=(iterator const &) = delete;
        iterator & operator=(iterator &&) = default;
        ~iterator() = default;

        constexpr reference operator*() const noexcept {
            return context{...};
        }

        constexpr pointer operator->() const noexcept {
            return context{...};
        }

        constexpr iterator & operator++() noexcept {
            assert(_host != nullptr);

            auto & parent = active_node();

            if (!parent.is_leaf()) {
                auto maybe_alt_child = parent.try_alt_child(); // simple strategy without checking the coverage?
                auto maybe_ref_child = parent.try_ref_child(); // could return an optional?
                if (maybe_ref_child.has_value()) {
                    parent = std::move(*maybe_ref_child);
                    if (maybe_alt_child.has_value())
                        visit(std::move(*maybe_alt_child));
                } else {
                    assert(maybe_alt_child.has_value())
                    parent = std::move(*alt_child);
                }
            } else { // go back to previous inner node.
                backtrack(); // to build an observable traverser, this must call some delegate function, or wraps it.
            }
            return *this;
        }

        constexpr void operator++(int) noexcept {
            ++(*this);
        }

    private:
        constexpr friend bool operator==(iterator const & rhs, std::default_sentinel_t) noexcept {
            return rhs._host->_active_branch.top().is_sink();
        }

        constexpr void visit(node_type && next_node) {
            _host->_active_branch.push_back(std::move(next_node));
        }

        constexpr node_type & active_node() noexcept {
            assert(_host != nullptr);
            return _host->_active_branch.top();
        }

        constexpr void backtrack() noexcept {
            assert(_host != nullptr);
            _host->_active_branch.pop();
        }
    };
}  // namespace libjst
