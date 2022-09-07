// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides serialisation concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/cereal.hpp>

#include <libcontrib/std/tag_invoke.hpp>

namespace libjst
{
    namespace _load
    {
        inline constexpr struct _cpo
        {
            template <typename object_t, typename archive_t>
            requires std::tag_invocable<_cpo, object_t, archive_t &>
            constexpr auto operator()(object_t &&obj, archive_t &archive) const
                -> std::tag_invoke_result_t<_cpo, object_t, archive_t &>
            {
                return std::tag_invoke(_cpo{}, (object_t &&) obj, archive);
            }
        } load;
    } // namespace _load
    using _load::load;

    namespace _save
    {
        inline constexpr struct _cpo
        {
            template <typename object_t, typename archive_t>
            requires std::tag_invocable<_cpo, object_t const &, archive_t &>
            constexpr auto operator()(object_t const &obj, archive_t &archive) const
                -> std::tag_invoke_result_t<_cpo, object_t const &, archive_t &>
            {
                return std::tag_invoke(_cpo{}, obj, archive);
            }
        } save;
    } // namespace _save
    using _save::save;

    namespace _load_extern
    {
        inline constexpr struct _cpo
        {
            template <typename archive_t, typename value_t>
            requires std::tag_invocable<_cpo, archive_t &, value_t const &>
            constexpr auto operator()(archive_t &archive, value_t const &value) const
                -> std::tag_invoke_result_t<_cpo, archive_t &, value_t const &>
            {
                return std::tag_invoke(_cpo{}, archive, value);
            }
        } load_extern;
    } // namespace _load_extern
    using _load_extern::load_extern;

    namespace _save_extern
    {
        inline constexpr struct _cpo
        {
            template <typename archive_t, typename value_t>
            requires std::tag_invocable<_cpo, archive_t &, value_t const &>
            constexpr auto operator()(archive_t &archive, value_t const &value) const
                -> std::tag_invoke_result_t<_cpo, archive_t &, value_t const &>
            {
                return std::tag_invoke(_cpo{}, archive, value);
            }
        } save_extern;
    } // namespace _save_extern
    using _save_extern::save_extern;
}  // namespace libjst
