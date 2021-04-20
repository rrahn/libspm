// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::search_state_manager_single.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <libjst/search/state_manager_concept.hpp>

namespace libjst
{

/*!\brief A search state manager handling a single state.
 * \implements libjst::search_state_manager
 *
 * \tparam search_state_t The type of the search state; must model std::semiregular.
 *
 * \details
 *
 * This manager is useful if the search only needs a single state during the execution, e.g. by searching in a
 * flat sequence.
 */
template <std::semiregular search_state_t>
class search_state_manager_single
{
private:
    search_state_t _state{}; //!\< The managed search state.
public:
    //!\brief The type of the managed search state.
    using state_type = search_state_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr search_state_manager_single() = default; //!< Default.
    constexpr search_state_manager_single(search_state_manager_single const &) = default; //!< Default.
    constexpr search_state_manager_single(search_state_manager_single &&) = default; //!< Default.
    constexpr search_state_manager_single & operator=(search_state_manager_single const &) = default; //!< Default.
    constexpr search_state_manager_single & operator=(search_state_manager_single &&) = default; //!< Default.
    ~search_state_manager_single() = default; //!< Default.
    //!\}

    //!\copydoc libjst::search_state_manager::state()
    constexpr state_type & state() noexcept
    {
        return _state;
    }

    //!\copydoc libjst::search_state_manager::state() const
    constexpr state_type const & state() const noexcept
    {
        return _state;
    }
};

}  // namespace libjst
