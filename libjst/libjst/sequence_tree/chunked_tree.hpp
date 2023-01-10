// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <cmath>

#include <libcontrib/closure_adaptor.hpp>

#include <libjst/sequence_tree/concept.hpp>
#include <libjst/sequence_tree/partial_tree.hpp>

namespace libjst
{
    template <typename wrapped_tree_t>
    class chunked_tree_impl : public std::ranges::view_base {
    private:

        using chunk_type = partial_tree<wrapped_tree_t const &>;
        using maybe_wrapped_t = std::conditional_t<std::is_lvalue_reference_v<wrapped_tree_t>,
                                                   std::reference_wrapper<std::remove_reference_t<wrapped_tree_t>>,
                                                   wrapped_tree_t>;

        class iterator;

        using reference = typename iterator::reference;
        using difference_type = typename iterator::difference_type;

        maybe_wrapped_t _wrappee;
        size_t _chunk_size{};

    public:

        constexpr chunked_tree_impl() = default;

        template <typename wrappee_t, std::unsigned_integral chunk_size_t>
            requires (!std::same_as<std::remove_cvref_t<wrappee_t>, chunked_tree_impl> &&
                      std::constructible_from<maybe_wrapped_t, wrappee_t>)
        explicit constexpr chunked_tree_impl(wrappee_t && wrappee, chunk_size_t chunk_size) noexcept :
            _wrappee{(wrappee_t &&)wrappee},
            _chunk_size{static_cast<size_t>(chunk_size)}
        {}

        chunked_tree_impl(chunked_tree_impl const &) = default;
        chunked_tree_impl(chunked_tree_impl &&) = default;
        chunked_tree_impl & operator=(chunked_tree_impl const &) = default;
        chunked_tree_impl & operator=(chunked_tree_impl &&) = default;
        ~chunked_tree_impl() = default;

        constexpr reference operator[](difference_type const step) const noexcept {
            return *(begin() + step);
        }

        constexpr iterator begin() const noexcept {
            return iterator{this, 0};
        }

        constexpr iterator end() const noexcept {
            return iterator{this, max_chunk_count()};
        }
    private:

        constexpr wrapped_tree_t const & base() const noexcept {
            if constexpr (std::is_lvalue_reference_v<wrapped_tree_t>)
                return _wrappee.get();
            else
                return _wrappee;
        }

        constexpr size_t max_chunk_count() const noexcept {
            return std::ceil(std::ranges::size(_wrappee.data().source())/_chunk_size);
        }
    };

    template <typename wrapped_tree_t>
    class chunked_tree_impl<wrapped_tree_t>::iterator {
    private:

        friend chunked_tree_impl;

        chunked_tree_impl const * _host{};
        size_t _chunk_idx{};

        constexpr iterator(chunked_tree_impl const * host, size_t initial_chunk_idx) noexcept :
            _host{host},
            _chunk_idx{initial_chunk_idx}
        {}

    public:

        using value_type = chunk_type;
        using reference = chunk_type;
        using pointer = void;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;

        iterator() = default;

        constexpr reference operator*() const noexcept {
            return chunk_type{_host->base(), chunk_offset(), _host->_chunk_size};
        }

        constexpr reference operator[](difference_type const step) const noexcept {
            return *(*this + step);
        }

        constexpr iterator & operator++() noexcept {
            ++_chunk_idx;
            return *this;
        }

        constexpr iterator operator++(int) noexcept {
            iterator tmp{*this};
            ++_chunk_idx;
            return tmp;
        }

        constexpr iterator & operator+=(difference_type const step) noexcept {
            _chunk_idx += step;
            return *this;
        }

        constexpr iterator & operator--() noexcept {
            --_chunk_idx;
            return *this;
        }

        constexpr iterator operator--(int) noexcept {
            iterator tmp{*this};
            --_chunk_idx;
            return tmp;
        }

        constexpr iterator & operator-=(difference_type const step) noexcept {
            _chunk_idx -= step;
            return *this;
        }

    private:

        constexpr auto chunk_offset() const noexcept {
            using offset_t = typename breakpoint::value_type;
            return static_cast<offset_t>(_chunk_idx * _host->_chunk_size);
        }

        constexpr friend iterator operator+(iterator lhs, difference_type const step) noexcept {
            return lhs += step;
        }

        constexpr friend iterator operator+(difference_type const step, iterator const & rhs) noexcept {
            return rhs + step;
        }

        constexpr friend iterator operator-(iterator lhs, difference_type const step) noexcept {
            return lhs -= step;
        }

        constexpr friend difference_type operator-(iterator const & lhs, iterator const & rhs) noexcept {
            return lhs._chunk_idx - rhs._chunk_idx;
        }

        constexpr friend bool operator==(iterator const &, iterator const &) noexcept = default;

        constexpr friend auto operator<=>(iterator const & lhs, iterator const & rhs) noexcept {
            return lhs._chunk_idx <=> rhs._chunk_idx;
        }
    };

    namespace _tree_adaptor {
        inline constexpr struct _chunk
        {
            template <typename wrapped_tree_t, std::unsigned_integral chunk_size_t>
            constexpr auto operator()(wrapped_tree_t && tree, chunk_size_t const chunk_size) const
                noexcept(std::is_nothrow_constructible_v<chunked_tree_impl<wrapped_tree_t>>)
                -> chunked_tree_impl<wrapped_tree_t>
            {
                using adapted_tree_t = chunked_tree_impl<wrapped_tree_t>;
                return adapted_tree_t{(wrapped_tree_t &&)tree, chunk_size};
            }

            template <std::unsigned_integral chunk_size_t>
            constexpr auto operator()(chunk_size_t const chunk_size) const
                noexcept(std::is_nothrow_invocable_v<std::tag_t<jst::contrib::make_closure>, chunk_size_t>)
                -> jst::contrib::closure_result_t<_chunk, chunk_size_t>
            { // we need to store the type that needs to be called later!
                return jst::contrib::make_closure(_chunk{}, chunk_size);
            }
        } chunk{};
    } // namespace _tree_adaptor

    using _tree_adaptor::chunk;
}  // namespace libjst
