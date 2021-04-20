// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::search_stack_observer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
/*!\interface libjst::search_stack_observer <>
 * \brief     The generic concept for stack observables.
 *
 * This concept describes the requirements a type must fulfil in order to notify it about stack events, i.e. push and
 * pop events.
 */
/*!\name Requirements for libjst::search_stack_observer
 * \brief You can expect these member functions on all types that model libjst::search_stack_observer.
 * \relates libjst::search_stack_observer
 * \{
 */
//!\cond
template <typename t>
concept search_stack_observer = requires (t & observer)
{
//!\endcond

    /*!\fn libjst::search_stack_observer::on_pop()
     * \brief Method to be invoked on a pop event.
     *
     * \attention This is a concept requirement, not an actual function (however types
     *            modeling this concept will provide an implementation).
     */
    { observer.on_pop() };

    /*!\fn libjst::search_stack_observer::on_push()
     * \brief Method to be invoked on a push event.
     *
     * \attention This is a concept requirement, not an actual function (however types
     *            modeling this concept will provide an implementation).
     */
    { observer.on_push() };
//!\cond
};
//!\endcond
//!\}
}  // namespace libjst
