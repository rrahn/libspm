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
#include <ranges>
#include <type_traits>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // CPO defintions
    // ----------------------------------------------------------------------------

    // position
    namespace _variant_position {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr auto operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
                -> std::tag_invoke_result_t<_cpo, variant_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
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
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr auto operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
                -> std::tag_invoke_result_t<_cpo, variant_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
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
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr auto operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
                -> std::tag_invoke_result_t<_cpo, variant_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
            }
        } deletion;
    }
    using _variant_deletion::deletion;

    template <typename variant_t>
    using variant_deletion_t = std::invoke_result_t<_variant_deletion::_cpo, variant_t>;

    // coverage
    namespace _variant_coverage {
        inline constexpr struct _cpo  {
            template <typename variant_t>
                requires std::tag_invocable<_cpo, variant_t const &>
            constexpr auto operator()(variant_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, variant_t const &>)
                -> std::tag_invoke_result_t<_cpo, variant_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
            }
        } coverage;
    }
    using _variant_coverage::coverage;

    template <typename variant_t>
    using variant_coverage_t = std::invoke_result_t<_variant_coverage::_cpo, variant_t>;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename type>
    concept sequence_variant = requires (type const & value)
    {
        { libjst::position(value) } -> std::integral;
        { libjst::deletion(value) } -> std::integral;
        { libjst::insertion(value) } -> std::ranges::forward_range;
    };

    template <typename type>
    concept covered_sequence_variant = sequence_variant<type> && requires (type const & value)
    {
        { libjst::coverage(value) } -> std::ranges::random_access_range;
    };

}  // namespace libjst
