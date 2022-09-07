// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides rooted rcs store implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libjst/variant/concept.hpp>

namespace libjst
{

    template <typename rcs_store_t>
        requires std::common_reference_with<
            std::ranges::range_value_t<typename rcs_store_t::variant_map_type> &,
            std::ranges::range_reference_t<typename rcs_store_t::variant_map_type>>
    class rooted_rcs_store
    {
    private:
        using base_variant_map_type = typename rcs_store_t::variant_map_type;

        class rooted_variant_map;

        std::reference_wrapper<rcs_store_t const> _wrappee;
        rooted_variant_map _rooted_variants{};

    public:

        using variant_map_type = rooted_variant_map;
        using source_type = typename rcs_store_t::source_type;

        rooted_rcs_store() = delete;
        explicit rooted_rcs_store(rcs_store_t const & wrappee) :
            _wrappee{std::cref(wrappee)},
            _rooted_variants{_wrappee.get().variants(), _wrappee.get().size()}
        {}

        constexpr rcs_store_t const & base() const noexcept {
            return _wrappee.get();
        }

        constexpr auto size() const noexcept -> decltype(_wrappee.get().size()){
            return _wrappee.get().size();
        }

        constexpr auto source() const noexcept -> decltype(_wrappee.get().source()){
            return _wrappee.get().source();
        }

        constexpr variant_map_type const & variants() const noexcept {
            return _rooted_variants;
        }

    };

    template <typename rcs_store_t>
        requires std::common_reference_with<
            std::ranges::range_value_t<typename rcs_store_t::variant_map_type> &,
            std::ranges::range_reference_t<typename rcs_store_t::variant_map_type>>
    class rooted_rcs_store<rcs_store_t>::rooted_variant_map {
    private:

        using value_type = std::ranges::range_value_t<base_variant_map_type>;

        template <bool>
        class iterator_type;

        value_type _root{};
        base_variant_map_type const & _wrappee{};

    public:

        using iterator = iterator_type<false>;
        using const_iterator = iterator_type<true>;

        rooted_variant_map(base_variant_map_type const & wrappee, size_t const coverage_size) : _wrappee{wrappee}
        {
            using coverage_t = libjst::variant_coverage_t<value_type>;

            coverage_t root_cov{};
            root_cov.resize(coverage_size, true);
            libjst::coverage(_root) = std::move(root_cov);
            // libjst::coverage(_root) = std::move(root_cov);
            // libjst::left_breakpoint(_root) = 0;
            // libjst::alt_sequence(_root) = libjst::alt_sequence_t<value_type>{};
            // libjst::breakpoint_span(_root) = 0;
        }

        iterator begin() noexcept {
            return iterator{std::ranges::begin(_wrappee), _root, -1};
        }

        const_iterator begin() const noexcept {
            return const_iterator{std::ranges::begin(_wrappee), _root, -1};
        }

        iterator end() noexcept {
            return iterator{std::ranges::end(_wrappee), _root, std::ranges::ssize(_wrappee)};
        }

        const_iterator end() const noexcept {
            return const_iterator{std::ranges::end(_wrappee), _root, std::ranges::ssize(_wrappee)};
        }
    };

    template <typename rcs_store_t>
        requires std::common_reference_with<
            std::ranges::range_value_t<typename rcs_store_t::variant_map_type> &,
            std::ranges::range_reference_t<typename rcs_store_t::variant_map_type>>
    template <bool is_const>
    class rooted_rcs_store<rcs_store_t>::rooted_variant_map::iterator_type {
    private:
        template <typename t>
        using maybe_const_t = std::conditional_t<is_const, t const, t>;

        using base_iterator = std::ranges::iterator_t<maybe_const_t<base_variant_map_type>>;

        friend rooted_variant_map;

    public:
        using value_type = std::iter_value_t<base_iterator>;
        using reference = std::common_reference_t<std::iter_reference_t<base_iterator>, maybe_const_t<value_type> &>;
        using difference_type = std::iter_difference_t<base_iterator>;
        using pointer_t = void;
        using iterator_category = std::iterator_traits<base_iterator>::iterator_category;

    private:

        base_iterator _it{};
        maybe_const_t<value_type> *_root{};
        difference_type _idx{};

        iterator_type(base_iterator it,
                      maybe_const_t<value_type> & root,
                      difference_type idx) : _it{std::move(it)}, _root{std::addressof(root)}, _idx{idx}
        {}

    public:

        iterator_type() = default;
        iterator_type(iterator_type<!is_const> other) requires is_const : _it{other._it}, _root{other.root}, _idx{other._idx}
        {}

        constexpr reference operator*() const noexcept {
            return is_root() ? static_cast<reference>(*_root) : *_it;
        }

        constexpr reference operator[](difference_type offset) const noexcept
            requires std::random_access_iterator<base_iterator>
        {
            iterator tmp = (*this + offset);
            return *tmp;
        }

        constexpr iterator_type & operator++() noexcept
            requires std::input_iterator<base_iterator>
        {
            if (_idx++ >= 0) ++_it;
            return *this;
        }

        constexpr iterator_type operator++(int) noexcept
            requires std::forward_iterator<base_iterator>
        {
            iterator_type tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator_type & operator--() noexcept
            requires std::bidirectional_iterator<base_iterator>
        {
            if (--_idx >= 0) --_it;
            return *this;
        }

        constexpr iterator_type operator--(int) noexcept
            requires std::bidirectional_iterator<base_iterator>
        {
            iterator_type tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator_type & operator+=(difference_type offset) noexcept
            requires std::random_access_iterator<base_iterator>
        {
            if (is_root() && offset > 0) {
                _idx += offset;
                _it += _idx;
            } else {
                _idx += offset;
                _it += offset;
            }
            return *this;
        }

        constexpr iterator_type & operator-=(difference_type offset) noexcept
            requires std::random_access_iterator<base_iterator>
        {
            _idx -= offset;
            if (_idx >= 0) {
                _it -= offset;
            } else {
                _it -= max(offset - 1, 0);
            }
        }

    private:

        constexpr bool is_root() const noexcept {
            return _idx == -1;
        }

        constexpr friend iterator_type operator+(iterator_type const & lhs, difference_type offset) noexcept
            requires std::random_access_iterator<base_iterator>
        {
            iterator_type tmp{lhs};
            return tmp += offset;
        }

        constexpr friend iterator_type operator+(difference_type offset, iterator_type const & rhs) noexcept
            requires std::random_access_iterator<base_iterator>
        {
            return rhs + offset;
        }

        constexpr friend iterator_type operator-(iterator_type lhs, difference_type offset) noexcept
            requires std::random_access_iterator<base_iterator>
        {
            iterator_type tmp{lhs};
            return tmp -= offset;
        }

        constexpr friend difference_type operator-(iterator_type const & lhs, iterator_type const & rhs) noexcept
            requires std::sized_sentinel_for<base_iterator, base_iterator>
        {
            return lhs._idx - rhs._idx;
        }

        constexpr friend bool operator==(iterator_type const & lhs, iterator_type const & rhs) noexcept {
            return lhs._idx == rhs._idx;
        }

        constexpr friend std::strong_ordering operator<=>(iterator_type const & lhs, iterator_type const & rhs) noexcept
        {
            return lhs._idx <=> rhs._idx;
        }
    };

}  // namespace libjst
