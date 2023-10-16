// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Implementation of tag_invoke utility to implement customisation point objects.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <utility>

namespace libjst
{
    // ----------------------------------------------------------------------------
    // Implementation of tag_invoke
    // ----------------------------------------------------------------------------

    /// @cond
    namespace _tag_invoke
    { // exposition only

        // Poison pill to prohibit ADL for tag_invoke further up in the namespace stack.
        void tag_invoke();

        template <typename tag_t, typename... args_t>
        concept tag_invocable = requires(tag_t &&tag, args_t &&...args) {
            tag_invoke(std::forward<tag_t>(tag), std::forward<args_t>(args)...);
        };

        template <typename tag_t, typename... args_t>
        concept nothrow_tag_invocable = tag_invocable<tag_t, args_t...> && requires(tag_t &&tag, args_t &&...args) {
            {
                tag_invoke(std::forward<tag_t>(tag), std::forward<args_t>(args)...)
            } noexcept;
        };

        template <typename tag_t, typename... args_t>
        inline constexpr bool is_nothrow_tag_invocable_v = false;

        template <typename tag_t, typename... args_t>
            requires nothrow_tag_invocable<tag_t, args_t...>
        inline constexpr bool is_nothrow_tag_invocable_v<tag_t, args_t...> = true;

        template <typename tag_t, typename... args_t>
        using tag_invoke_result_t = decltype(tag_invoke(std::declval<tag_t>(), std::declval<args_t>()...));

        template <typename tag_t, typename... args_t>
        struct tag_invoke_result
        {
            using type = tag_invoke_result_t<tag_t, args_t...>; // present if and only if tag_invocable<tag_t, args_t...> is true
        };

        template <auto &tag_v>
        using tag_t = std::remove_cvref_t<decltype(tag_v)>;

        struct fn
        {
            template <typename tag_t, typename... args_t>
            constexpr auto operator()(tag_t &&tag, args_t &&...args) const
                noexcept(noexcept(tag_invoke(std::forward<tag_t>(tag), std::forward<args_t>(args)...)))
                    -> decltype(tag_invoke(std::forward<tag_t>(tag), std::forward<args_t>(args)...))
            {
                return tag_invoke(std::forward<tag_t>(tag), std::forward<args_t>(args)...);
            }
        };

        inline namespace _cpo
        {
            inline constexpr _tag_invoke::fn tag_invoke{};
        } // namespace _cpo

    } // namespace _tag_invoke
    /// @endcond

    using namespace _tag_invoke::_cpo;
    using _tag_invoke::is_nothrow_tag_invocable_v;
    using _tag_invoke::nothrow_tag_invocable;
    using _tag_invoke::tag_invocable;
    using _tag_invoke::tag_invoke_result;
    using _tag_invoke::tag_invoke_result_t;
    using _tag_invoke::tag_t;

} // namespace libjst
