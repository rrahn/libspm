// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides base implementation of a reversed referentially compressed multisequence (rcms).
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <libcontrib/type_traits.hpp>
#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/coverage/concept.hpp>
#include <libjst/variant/concept.hpp>
#include <libjst/variant/breakpoint_reversed.hpp>

namespace libjst
{
    // requires is_object_v<source_t>
    // requires breakpoint_coverage<coverage_t>
    template <typename rcms_t>
    class compressed_multisequence_reversed {
    private:

        using wrapped_source_t = decltype(std::declval<rcms_t const &>().source());
        using wrapped_iterator = std::ranges::iterator_t<rcms_t const &>;
        using coverage_type = libjst::variant_coverage_t<std::iter_reference_t<wrapped_iterator>>;
        using coverage_domain_type = libjst::coverage_domain_t<coverage_type>;

        class iterator_impl;
        class delta_proxy;

        rcms_t const & _wrappee;

    public:

        using source_type = std::ranges::reverse_view<std::views::all_t<wrapped_source_t>>;
        using iterator = iterator_impl;
        using value_type = std::iter_value_t<iterator>;

        compressed_multisequence_reversed() = default;
        explicit compressed_multisequence_reversed(rcms_t const & wrappee) :
            _wrappee{wrappee}
        {
        }

        constexpr source_type source() const noexcept {
            return _wrappee.source() | std::views::all | std::views::reverse;
        }

        constexpr size_t size() const noexcept {
            return _wrappee.size();
        }

        constexpr coverage_domain_type const & coverage_domain() const noexcept {
            return _wrappee.coverage_domain();
        }

        constexpr iterator begin() const noexcept {
            return iterator{_wrappee.end(), std::ranges::size(source())};
        }

        constexpr iterator end() const noexcept {
            return iterator{_wrappee.begin(), std::ranges::size(source())};
        }

    private:

    };

    template <typename rcms_t>
    class compressed_multisequence_reversed<rcms_t>::iterator_impl {

        friend compressed_multisequence_reversed;

        wrapped_iterator _reverse_it{};
        size_t _source_size{};

        explicit constexpr iterator_impl(wrapped_iterator reverse_it, size_t source_size) noexcept :
            _reverse_it{std::move(reverse_it)},
            _source_size{source_size}
        {}

    public:

        using value_type = delta_proxy;
        using reference = value_type;
        using difference_type = std::iter_difference_t<wrapped_iterator>;
        using pointer = void;
        using iterator_category = std::random_access_iterator_tag;

        constexpr iterator_impl() = default;

        constexpr reference operator*() const noexcept {
            return reference{std::ranges::prev(_reverse_it), _source_size};
        }

        constexpr reference operator[](difference_type const step) const noexcept {
            return *(*this + step);
        }

        constexpr iterator_impl & operator++() noexcept {
            --_reverse_it;
            return *this;
        }

        constexpr iterator_impl operator++(int) noexcept {
            iterator_impl tmp{*this};
            ++(*this);
            return tmp;
        }

        constexpr iterator_impl & operator+=(difference_type const step) noexcept {
            _reverse_it -= step;
            return *this;
        }

        constexpr iterator_impl & operator--() noexcept {
            ++_reverse_it;
            return *this;
        }

        constexpr iterator_impl operator--(int) noexcept {
            iterator_impl tmp{*this};
            --(*this);
            return tmp;
        }

        constexpr iterator_impl & operator-=(difference_type const step) noexcept {
            _reverse_it += step;
            return *this;
        }

    private:

        constexpr friend iterator_impl operator+(iterator_impl lhs, difference_type const step) noexcept {
            return lhs += step;
        }

        constexpr friend iterator_impl operator+(difference_type const step, iterator_impl rhs) noexcept {
            return rhs + step;
        }

        constexpr friend iterator_impl operator-(iterator_impl lhs, difference_type const step) noexcept {
            return lhs -= step;
        }

        constexpr friend difference_type operator-(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return rhs._reverse_it - lhs._reverse_it;
        }

        constexpr friend bool operator==(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._reverse_it == rhs._reverse_it;
        }

        constexpr friend std::strong_ordering operator<=>(iterator_impl const & lhs, iterator_impl const & rhs) noexcept {
            return lhs._reverse_it <=> rhs._reverse_it;
        }
    };

    template <typename rcms_t>
    class compressed_multisequence_reversed<rcms_t>::delta_proxy {

        friend compressed_multisequence_reversed;

        using coverage_t = libjst::variant_coverage_t<std::iter_reference_t<wrapped_iterator>>;
        using position_t = libjst::variant_position_t<std::iter_reference_t<wrapped_iterator>>;

        wrapped_iterator _reverse_it{};
        size_t _source_size{};

        explicit constexpr delta_proxy(wrapped_iterator reverse_it,
                                       size_t source_size) noexcept :
            _reverse_it{std::move(reverse_it)},
            _source_size{source_size}
        {}

    public:

        delta_proxy() = delete;

        constexpr breakpoint_end get_breakpoint_end() const noexcept {
            breakpoint_end tmp = (*_reverse_it).get_breakpoint_end();
            return (tmp == breakpoint_end::low) ? breakpoint_end::high : breakpoint_end::low;
        }

        constexpr std::optional<wrapped_iterator> jump_to_mate() const noexcept {
            return (*_reverse_it).jump_to_mate();
        }

    private:

        friend constexpr position_t tag_invoke(std::tag_t<libjst::low_breakend>, delta_proxy me) noexcept
        {
            return libjst::low_breakend(libjst::get_breakpoint(me));
        }

        friend constexpr position_t tag_invoke(std::tag_t<libjst::high_breakend>, delta_proxy me) noexcept
        {
            return libjst::high_breakend(libjst::get_breakpoint(me));
        }

        friend constexpr breakpoint_reversed tag_invoke(std::tag_t<libjst::get_breakpoint>, delta_proxy me) noexcept
        {
            return breakpoint_reversed{libjst::get_breakpoint(*me._reverse_it), me._source_size};
        }

        friend constexpr auto tag_invoke(std::tag_t<libjst::alt_sequence>, delta_proxy me) noexcept
        {
            return libjst::alt_sequence(*me._reverse_it) | std::views::reverse;
        }

        friend constexpr position_t tag_invoke(std::tag_t<libjst::position>, delta_proxy me) noexcept
        {
            return me._source_size - libjst::position(*me._reverse_it);
        }

        template <typename cpo_t>
            requires std::tag_invocable<cpo_t, std::iter_reference_t<wrapped_iterator>>
        friend constexpr auto tag_invoke(cpo_t cpo, delta_proxy me)
            noexcept(std::is_nothrow_tag_invocable_v<cpo_t, std::iter_reference_t<wrapped_iterator>>)
            -> std::tag_invoke_result_t<cpo_t, std::iter_reference_t<wrapped_iterator>>
        {
            return cpo(*me._reverse_it);
        }
    };

}  // namespace libjst
