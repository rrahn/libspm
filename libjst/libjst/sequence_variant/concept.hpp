// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for sequence variants.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <ranges>
#include <type_traits>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Operation CPOs for variants
    // ----------------------------------------------------------------------------

    // position
    namespace _variant_position {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t && var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> std::tag_invoke_result_t<_cpo, variant_t>
            {
                return std::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } position;
    } // namespace _variant_position
    using _variant_position::position;

    template <typename variant_t>
    using variant_position_t = std::invoke_result_t<_variant_position::_cpo, variant_t>;

    // insertion
    namespace _variant_insertion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> std::tag_invoke_result_t<_cpo, variant_t>
            {
                return std::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } insertion;
    } // namespace insertion
    using _variant_insertion::insertion;

    template <typename variant_t>
    using variant_insertion_t = std::invoke_result_t<_variant_insertion::_cpo, variant_t>;

    // deletion
    namespace _variant_deletion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> std::tag_invoke_result_t<_cpo, variant_t>
            {
                return std::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } deletion;
    }
    using _variant_deletion::deletion;

    template <typename variant_t>
    using variant_deletion_t = std::invoke_result_t<_variant_deletion::_cpo, variant_t>;

    //is_deletion
    namespace _is_deletion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return std::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!std::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (std::tag_invocable<std::tag_t<libjst::insertion>, variant_t const &> &&
                              std::tag_invocable<std::tag_t<libjst::deletion>, variant_t const &>)
                    return libjst::deletion(var) > 0 && std::ranges::size(libjst::insertion(var)) == 0;
                else
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_deletion(variant) found!");
            }

        } is_deletion;
    }
    using _is_deletion::is_deletion;

    //is_insertion
    namespace _is_insertion {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return std::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!std::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (std::tag_invocable<std::tag_t<libjst::insertion>, variant_t const &> &&
                              std::tag_invocable<std::tag_t<libjst::deletion>, variant_t const &>)
                    return libjst::deletion(var) == 0 && std::ranges::size(libjst::insertion(var)) > 0;
                else
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_insertion(variant) found!");
            }

        } is_insertion;
    }
    using _is_insertion::is_insertion;

    //is_replacement
    namespace _is_replacement {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr bool operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
            {
                return std::tag_invoke(_cpo{}, var);
            }

            template <typename variant_t>
                requires (!std::tag_invocable<_cpo, variant_t const &>)
            constexpr bool operator()(variant_t const & var) const noexcept
            {
                if constexpr (std::tag_invocable<std::tag_t<libjst::insertion>, variant_t const &> &&
                              std::tag_invocable<std::tag_t<libjst::deletion>, variant_t const &>)
                {
                    return libjst::deletion(var) == std::ranges::size(libjst::insertion(var)) &&
                           libjst::deletion(var) > 0;
                } else {
                    static_assert(std::same_as<variant_t, void>, "No valid default for is_replacement(variant) found!");
                }
            }

        } is_replacement;
    }
    using _is_replacement::is_replacement;

    // coverage
    namespace _variant_coverage {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t>
            constexpr auto operator()(variant_t &&var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t>)
                -> std::tag_invoke_result_t<_cpo, variant_t>
            {
                return std::tag_invoke(_cpo{}, (variant_t &&)var);
            }
        } coverage;
    }
    using _variant_coverage::coverage;

    template <typename variant_t>
    using variant_coverage_t = std::remove_cvref_t<std::invoke_result_t<_variant_coverage::_cpo, variant_t>>;

    // ----------------------------------------------------------------------------
    // Operation CPOs for variant stores
    // ----------------------------------------------------------------------------

    namespace _insert {
        inline constexpr struct _cpo  {
            template <typename variant_store_t, typename variant_t>
                requires std::tag_invocable<_cpo, variant_store_t &, variant_t>
            constexpr auto operator()(variant_store_t & store, variant_t && variant) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_store_t &, variant_t>)
                -> std::tag_invoke_result_t<_cpo, variant_store_t &, variant_t>
            {
                return std::tag_invoke(_cpo{}, store, (variant_t &&)variant);
            }
        } insert;
    }
    using _insert::insert;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename variant_t>
    concept sequence_variant = requires
    (variant_t const & variant)
    {
        { libjst::position(variant) } -> std::integral;
        { libjst::deletion(variant) } -> std::integral;
        { libjst::insertion(variant) } -> std::ranges::forward_range;
        { libjst::is_replacement(variant) } -> std::same_as<bool>;
        { libjst::is_insertion(variant) } -> std::same_as<bool>;
        { libjst::is_deletion(variant) } -> std::same_as<bool>;
    };

    template <typename variant_t>
    concept covered_sequence_variant = sequence_variant<variant_t> && requires
    (variant_t const & variant)
    {
        { libjst::coverage(variant) } -> std::ranges::random_access_range;
    };

    template <typename store_t>
    concept sequence_variant_store = std::ranges::random_access_range<store_t>
    && requires
    (store_t const & store, std::ranges::range_value_t<store_t> const & variant)
    {
        requires sequence_variant<std::ranges::range_value_t<store_t>>;
        requires sequence_variant<std::ranges::range_reference_t<store_t>>;

        // { store.insert(variant) } -> std::random_access_iterator;
    }
    ;

    template <typename store_t>
    concept covered_sequence_variant_store = sequence_variant_store<store_t>
    && requires
    {
        requires covered_sequence_variant<std::ranges::range_value_t<store_t>>;
        requires covered_sequence_variant<std::ranges::range_reference_t<store_t>>;
    }
    ;

}  // namespace libjst
