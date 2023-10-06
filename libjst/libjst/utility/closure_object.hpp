// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides closure adaptor for all kind of objects.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

namespace libjst
{

    namespace detail {

        struct closure_base{};

        struct nested_closure
        {
            template <typename target_t, typename inner_closure_t, typename outer_closure_t>
            constexpr auto operator()(target_t && target, inner_closure_t && inner, outer_closure_t && outer) const
                noexcept(std::is_nothrow_invocable_v<inner_closure_t, target_t> &&
                         std::is_nothrow_invocable_v<outer_closure_t, std::invoke_result_t<inner_closure_t, target_t>>)
                -> std::invoke_result_t<outer_closure_t, std::invoke_result_t<inner_closure_t, target_t>>
            {
                return std::invoke((outer_closure_t &&)outer,
                                   std::invoke((inner_closure_t &&)inner, (target_t &&) target));
            }
        };

        template <typename enclosee_t, typename ...args_t>
        class closure : closure_base
        {
        private:
            [[no_unique_address]] enclosee_t _enclosee; // the target object of the closure operation
            std::tuple<args_t...> _enclosed_args;

        public:

            template <typename _enclosee_t, typename ..._args_t>
                requires (!std::same_as<std::remove_cvref_t<_enclosee_t>, closure>)
            constexpr explicit closure(_enclosee_t &&enclosee, _args_t &&...args)
                noexcept(std::is_nothrow_constructible_v<enclosee_t, _enclosee_t> &&
                         std::is_nothrow_constructible_v<std::tuple<args_t...>, _args_t...>) :
                _enclosee{(_enclosee_t &&)enclosee},
                _enclosed_args{(_args_t &&)args...}
            {}

            // Copy arguments if this is an lvalue!
            template <typename target_t>
                requires (!std::is_base_of_v<closure_base, std::remove_cvref_t<target_t>> &&
                          std::invocable<enclosee_t const &, target_t, args_t const &...>)
            constexpr auto operator()(target_t && target) const &
                noexcept(std::is_nothrow_invocable_v<enclosee_t const &, target_t, args_t const &...>)
                -> std::invoke_result_t<enclosee_t const &, target_t, args_t const &...>
            { // args are copied from result
                auto copy_apply = [&] (auto const &...enclosed_args) {
                    return std::invoke(_enclosee, (target_t &&)target, enclosed_args...);
                };
                return std::apply(copy_apply, _enclosed_args); // we need to apply now the arguments
            }

            // Move arguments if it is a temporary!
            template <typename target_t>
                requires (!std::is_base_of_v<closure_base, std::remove_cvref_t<target_t>> &&
                          std::invocable<enclosee_t, target_t, args_t...>)
            constexpr auto operator()(target_t && target) &&
                noexcept(std::is_nothrow_invocable_v<enclosee_t, target_t, args_t...>)
                -> std::invoke_result_t<enclosee_t, target_t, args_t...>
            {
                auto move_apply = [&] (auto &&...enclosed_args) {
                    return std::invoke((enclosee_t &&)_enclosee ,(target_t &&)target, std::move(enclosed_args)...);
                };
                return std::apply(move_apply, std::move(_enclosed_args)); // we need to apply now the arguments
            }

            // Apply the enclosed arguments with the target on the enclosee.
            template <typename target_t, typename closure_t>
                requires (!std::is_base_of_v<closure_base, std::remove_cvref_t<target_t>> &&
                          std::same_as<std::remove_cvref_t<closure_t>, closure> &&
                          std::invocable<closure_t, target_t>)
            friend constexpr auto operator|(target_t && target, closure_t && me)
                noexcept(std::is_nothrow_invocable_v<closure_t, target_t>)
                -> std::invoke_result_t<closure_t, target_t>
            {
                return std::invoke((closure_t &&)me, (target_t &&)target);
            }

            // Produce a nested closure object calling the rhs closure on the lhs closure after invoking it with the
            // target.
            template <typename fst_closure_t, typename snd_closure_t>
                requires (std::is_base_of_v<closure_base, std::remove_cvref_t<fst_closure_t>> &&
                          std::same_as<std::remove_cvref_t<snd_closure_t>, closure>)
            friend constexpr auto operator|(fst_closure_t && fst_closure, snd_closure_t && me)
                noexcept(std::is_nothrow_constructible_v<closure<nested_closure,
                                                                 std::decay_t<fst_closure_t>,
                                                                 std::decay_t<snd_closure_t>>>)
                -> closure<nested_closure, std::decay_t<fst_closure_t>, std::decay_t<snd_closure_t>>
            {
                // return a new closure object that resolves the nested closures in the correct order.
                using closure_t = closure<nested_closure, std::decay_t<fst_closure_t>, std::decay_t<snd_closure_t>>;
                return closure_t{nested_closure{}, (fst_closure_t &&)fst_closure, (snd_closure_t &&)me};
            }

        };
    } // namespace detail

    namespace _make_closure {
        inline constexpr struct _fn
        {
            template <typename enclosee_t, typename ...args_t>
            using closure_t = detail::closure<enclosee_t, std::decay_t<args_t>...>;

            // now we are creating the closure from a valid enclosee object
            template <typename enclosee_t, typename ...args_t>
            constexpr auto operator()(enclosee_t && enclosee, args_t &&...args) const
                noexcept(std::is_nothrow_constructible_v<closure_t<enclosee_t, args_t...>, enclosee_t, args_t...>)
                -> closure_t<enclosee_t, args_t...>
            {
                return closure_t<enclosee_t, args_t...>{(enclosee_t &&)enclosee, (args_t &&)args...};
            }
        } make_closure;
    } // namespace _closure

    using _make_closure::make_closure;

    template <typename enclosee_t, typename ...args_t>
    using closure_result_t = std::invoke_result_t<_make_closure::_fn, enclosee_t, args_t...>;

}  // namespace libjst
