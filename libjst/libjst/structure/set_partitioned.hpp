// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the partitioned journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/structure/concept_jst.hpp>
#include <libjst/structure/concept_serialiser.hpp>

namespace libjst
{
    namespace detail {

        template <traversable_journaled_sequence_tree jst_t>
        class jst_partitioned : private traversable_jst_base
        {
        private:

            [[no_unique_address]] jst_t const *_jst{};
            size_t _bin_index{};

        public:
            jst_partitioned() = default;
            explicit jst_partitioned(jst_t const * jst, size_t const bin_index) :
                _jst{jst},
                _bin_index{bin_index}
            {}

        private:

            // delegate to all getter interfaces of the wrapped component
            template <typename cpo_t>
                requires std::invocable<cpo_t, jst_t const &>
            constexpr friend auto tag_invoke(cpo_t cpo, jst_partitioned const &me)
                noexcept(std::is_nothrow_invocable_v<cpo_t, jst_t const &>)
                -> std::invoke_result_t<cpo_t, jst_t const &>
            {
                assert(me._jst != nullptr);
                return cpo(*me._jst);
            }
        };
    } // namespace detail

    // wrapper interface for the partial jst:
    template <traversable_journed_sequence_tree jst_t>
    class journaled_sequence_tree_partitioned : private traversable_jst_base
    {
    private:
        using store_t = variant_store_t<jst_t const &>;

        [[no_unique_address]] jst_t const *_jst{};
        size_t _bin_count{};

        // what is he returning and what kind of node type to we get when we access the node?
        class iterator;

    public:

        journaled_sequence_tree_partitioned() = default;
        explicit journaled_sequence_tree_partitioned(jst_t const & jst, size_t const bin_count) :
            _jst{std::addressof(jst)},
            _bin_count{bin_count}
        {}

        iterator begin() const noexcept
        {

        }

        iterator end() const noexcept
        {

        }
    private:
        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load>,
                                         journaled_sequence_tree_partitioned & me,
                                         archive_t & archive)
        {
            libjst::load_extern(archive, *me._jst);
            archive(me._bin_count);
        }

        template <typename archive_t>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save>,
                                         journaled_sequence_tree_partitioned const & me,
                                         archive_t & archive)
        {
            libjst::save_extern(archive, *me._jst);
            archive(me._bin_count);
        }
    };

    template <traversable_journed_sequence_tree _jst>
    class journaled_sequence_tree_partitioned<jst_t>::iterator
    {
    private:
        jst_t const * _jst{};
        std::ptrdiff_t _bin_idx{};

    public:

        using value_type = detail::jst_partitioned<jst_t>;
        using reference = value_type;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using iterator_category = std::random_access_iterator_tag;

        iterator() = default;
        iterator(_jst const * jst, std::ptrdiff_t bin_idx) noexcept : _jst{jst}, _bin_idx{bin_idx}
        {
        }

        // proxy type!
        reference operator*() const noexcept
        {
            return reference{_jst, _bin_idx};
        }

        iterator & operator++() noexcept
        {
            ++_bin_idx;
            return *this;
        }

        iterator operator++(int) noexcept
        {
            iterator tmp{*this};
            ++_bin_idx;
            return tmp;
        }

        iterator & operator+=(difference_type const offset) noexcept
        {
            _bin_idx += offset;
            return *this;
        }

        iterator operator+(difference_type const offset) const noexcept
        {
            iterator tmp{*this};
            tmp += offset;
            return tmp;
        }

        friend iterator operator+(difference_type const offset, iterator const & rhs) noexcept
        {
            return rhs + offset;
        }

        iterator & operator--() noexcept
        {
            --_bin_idx;
            return *this;
        }

        iterator operator--(int) noexcept
        {
            iterator tmp{*this};
            --_bin_idx;
            return tmp;
        }

        iterator & operator-=(difference_type const offset) noexcept
        {
            _bin_idx -= offset;
            return *this;
        }

        iterator operator-(difference_type const offset) const noexcept
        {
            iterator tmp{*this};
            tmp -= offset;
            return tmp;
        }

        difference_type operator-(iterator const & rhs) const noexcept
        {
            return _bin_idx - rhs._bin_idx;
        }

        constexpr bool operator==(iterator const & rhs) const noexcept =  default;

        constexpr std::strong_ordering operator<=>(iterator const & rhs) const noexcept
        {
            return _bin_idx <=> rhs._bin_idx;
        }
    };
}  // namespace libjst
