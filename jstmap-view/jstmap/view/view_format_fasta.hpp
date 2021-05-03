// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the functionality to view a particular haplotype of the jst.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <jstmap/view/global_types.hpp>
#include <jstmap/view/options.hpp>

namespace jstmap
{

//!\brief Displays haplotype for the specified index in fasta format.
void view_as_format(jst_t const & jst, size_t const index);

} // namespace jstmap
