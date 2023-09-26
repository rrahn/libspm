// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides composite class to wrap multiple invocable function objects.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{

    /**
     * @brief A type that can conveniently inherit multiple invocables and acts as a union over them.
     * @tparam invocable_ts The inherited types whose function call-operator gets invoked.
     */
    template <typename... invocable_ts>
    struct multi_invocable : invocable_ts...
    {
        /// @brief Inherited function call operators.
        using invocable_ts::operator()...;
    };

    /** @name Deduction guides */
    ///@{
    /**
     * @relates multi_invocable
     *
     * @tparam invocable_ts The inherited types whose function call-operator gets invoked.
     *
     * Deduces the template parameters from the given constructor arguments.
     */
    template <typename... invocable_ts>
    multi_invocable(invocable_ts...) -> multi_invocable<invocable_ts...>;
    ///@}

}  // namespace libjst
