// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a wrapper to encode indel variants in an rcms object.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <cereal/types/variant.hpp>

namespace libjst
{

    template <typename mate_iterator_t>
    class deletion_element {
        mate_iterator_t _mate_it{};
    public:

        using value_type = mate_iterator_t;

        constexpr deletion_element() = default;
        constexpr deletion_element(value_type mate_it)
            noexcept(std::is_nothrow_move_constructible_v<value_type>) : _mate_it{mate_it}
        {}

        constexpr value_type const & value() const noexcept {
            return _mate_it;
        }

        constexpr value_type & value() noexcept {
            return _mate_it;
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <typename archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(value());
        }

        template <typename archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(value());
        }

    };

    template <typename source_t>
    class insertion_element {
        source_t _insertion{};

    public:

        using value_type = source_t;

        constexpr insertion_element() = default;
        constexpr insertion_element(value_type insertion)
            noexcept(std::is_nothrow_move_constructible_v<value_type>) : _insertion{std::move(insertion)}
        {}

        constexpr value_type const & value() const noexcept {
            return _insertion;
        }

        constexpr value_type & value() noexcept {
            return _insertion;
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <typename archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(value());
        }

        template <typename archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(value());
        }

    };

    template <typename deletion_t, typename insertion_t>
    class indel_variant {

        using indel_type = std::variant<deletion_t, insertion_t>;

        indel_type _indel{};

    public:

        using value_type = indel_type;

        constexpr indel_variant() = default;
        constexpr indel_variant(value_type indel) noexcept(std::is_nothrow_move_constructible_v<value_type>) :
            _indel{std::move(indel)}
        {}

        constexpr value_type const & value() const noexcept {
            return _indel;
        }

        constexpr value_type & value() noexcept {
            return _indel;
        }

        template <typename visitor_t>
        constexpr auto visit(visitor_t && visitor) const
            noexcept(noexcept(std::visit((visitor_t&&)visitor, std::declval<indel_type const &>())))
            -> decltype(std::visit((visitor_t&&)visitor, std::declval<indel_type const &>()))
        {
            return std::visit((visitor_t&&)visitor, _indel);
        }

        // ----------------------------------------------------------------------------
        // Serialisation
        // ----------------------------------------------------------------------------

        template <typename archive_t>
        void load(archive_t & iarchive)
        {
            iarchive(_indel);
        }

        template <typename archive_t>
        void save(archive_t & oarchive) const
        {
            oarchive(_indel);
        }

    };
}  // namespace libjst
