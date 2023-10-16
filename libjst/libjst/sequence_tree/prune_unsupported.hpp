// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides prune tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/sequence_tree/prune_tree.hpp>

namespace libjst
{
    // namespace _tree_adaptor {
    //     inline constexpr struct _prune_unsupported
    //     {
    //         template <typename covered_tree_t, typename ...args_t>
    //         constexpr auto operator()(covered_tree_t && tree, args_t &&... args) const
    //             noexcept(std::is_nothrow_constructible_v<
    //                         prune_tree_impl<std::remove_reference_t<covered_tree_t>>, args_t...>)
    //             -> prune_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>
    //         {
    //             using adapted_tree_t = prune_tree_impl<std::remove_reference_t<covered_tree_t>, args_t...>;
    //             return adapted_tree_t{(covered_tree_t &&)tree, (args_t &&)args...};
    //         }

    //         template <typename ...args_t>
    //         constexpr auto operator()(args_t &&... args) const
    //             noexcept(std::is_nothrow_invocable_v<libjst::tag_t<libjst::make_closure>, args_t...>)
    //             -> libjst::closure_result_t<_prune_unsupported, args_t...>
    //         { // we need to store the type that needs to be called later!
    //             return libjst::make_closure(_prune_unsupported{}, (args_t &&)args...);
    //         }
    //     } prune_unsupported{};
    // } // namespace _tree_adaptor

    inline constexpr auto prune_unsupported = _tree_adaptor::prune;
}  // namespace libjst
