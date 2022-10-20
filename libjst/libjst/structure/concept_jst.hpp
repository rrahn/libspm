// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for journaled sequence tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

#include <libjst/variant/concept.hpp>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Operation CPOs
    // ----------------------------------------------------------------------------

    // base_sequence
    namespace _base_sequence {
        inline constexpr struct _cpo  {
            template <typename jst_t>
                requires std::tag_invocable<_cpo, jst_t const &>
            constexpr auto operator()(jst_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, jst_t const &>)
                -> std::tag_invoke_result_t<_cpo, jst_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
            }
        } base_sequence;
    } // namespace _base_sequence
    using _base_sequence::base_sequence;

    template <typename jst_t>
    using base_sequence_t = std::remove_cvref_t<std::invoke_result_t<_base_sequence::_cpo, jst_t>>;

    // variant_store
    namespace _variant_store {
        inline constexpr struct _cpo  {
            template <typename jst_t>
                requires std::tag_invocable<_cpo, jst_t const &>
            constexpr auto operator()(jst_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, jst_t const &>)
                -> std::tag_invoke_result_t<_cpo, jst_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
            }
        } variant_store;
    } // namespace _variant_store
    using _variant_store::variant_store;

    template <typename jst_t>
    using variant_store_t = std::remove_cvref_t<std::invoke_result_t<_variant_store::_cpo, jst_t>>;

    // size
    namespace _size {
        inline constexpr struct _cpo  {
            template <typename jst_t>
                requires std::tag_invocable<_cpo, jst_t const &>
            constexpr auto operator()(jst_t const & var) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, jst_t const &>)
                -> std::tag_invoke_result_t<_cpo, jst_t const &>
            {
                return std::tag_invoke(_cpo{}, var);
            }
        } size;
    } // namespace _size
    using _size::size;

    // path
    namespace _path {
        inline constexpr struct _cpo  {
            template <typename jst_t>
                requires std::tag_invocable<_cpo, jst_t const &>
            constexpr auto operator()(jst_t const & jst) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, jst_t const &>)
                -> std::tag_invoke_result_t<_cpo, jst_t const &>
            {
                return std::tag_invoke(_cpo{}, jst);
            }
        } path;
    } // namespace _path
    using _path::path;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------
    inline namespace current {
    template <typename jst_t>
    concept journaled_sequence_tree_c = requires (jst_t const & jst)
    {
        typename base_sequence_t<jst_t>;
        typename variant_store_t<jst_t>;
        // TODO: common type between reference and variant sequence

        { libjst::variant_store(jst) } -> libjst::covered_sequence_variant_store;
        { libjst::base_sequence(jst) } -> std::ranges::random_access_range;
        { libjst::size(jst) } -> std::unsigned_integral;
    };

    // Being traversable is a semantic property not testable by a concept definition.
    // Correspondingly, we use a tag to denote jst types that are traversable, i.e. only jst classes
    // that are derived from `traversable_jst_base` are traversable.
    struct traversable_jst_base{};

    template <typename jst_t>
    concept traversable_journaled_sequence_tree = journaled_sequence_tree_c<jst_t> &&
                                                  std::is_base_of_v<traversable_jst_base, std::remove_cvref_t<jst_t>> &&
    requires (jst_t && jst)
    {
        libjst::path((jst_t &&)jst);
    };
    } // namespace current
}  // namespace libjst
