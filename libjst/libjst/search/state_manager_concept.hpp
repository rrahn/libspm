// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::search_state_manager.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
/*!\interface libjst::search_state_manager <>
 * \brief The generic concept for state manager.
 *
 * This concept describes the requirements a type must fulfil in order to manage the state for a search algorithm.
 */
/*!\name Requirements for libjst::search_state_manager
 * \brief You can expect these member functions on all types that model libjst::search_state_manager.
 * \relates libjst::search_state_manager
 * \{
 */
//!\cond
template <typename t>
concept search_state_manager = requires (t & manager, t const & const_manager)
{
//!\endcond
    /*!\typedef state_type
     * \brief The type of state managed by the object.
     */
    typename t::state_type;

    /*!\fn state_reference libjst::search_state_manager::state()
     * \brief Returns a reference to the currently managed state.
     *
     * \details
     *
     * The returned type must model `std::assignable_from<state_reference, state_type>`.
     *
     * \attention This is a concept requirement, not an actual function (however types
     *            modeling this concept will provide an implementation).
     */
    { manager.state() } -> std::assignable_from<typename t::state_type>;

    /*!\fn state_const_reference libjst::search_state_manager::state() const
     * \brief Returns a const reference to the currently managed state.
     *
     * \details
     *
     * The returned type must model `std::convertible_to<state_const_reference, state_type>`.
     *
     * \attention This is a concept requirement, not an actual function (however types
     *            modeling this concept will provide an implementation).
     */
    { const_manager.state() } -> std::convertible_to<typename t::state_type>;
//!\cond
};
//!\endcond
//!\}

}  // namespace libjst
