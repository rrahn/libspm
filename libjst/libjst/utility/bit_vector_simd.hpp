// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/utility/bit_vector.hpp>

namespace libjst
{

/*!\brief A vectorised version of the bit vector.
 *
 * \tparam allocator_t The type of the allocator to use. Defaults to seqan3::aligned_allocator.
 *
 * \details
 *
 * Implements a bit vector with vectorised operations binary operations. The host memory should be aligned to a
 * multiple of the simd vector alignment specification. It always reserves more memory in order to read a full
 * vector from the underlying host memory without accessing memory that was not reserved by the bit vector.
 *
 * The reference type is a special proxy that provides access to a single bit. Note that it is not a real reference
 * but can be converter to a bool or assigned from a bool.
 */
template <typename allocator_t = std::allocator<uint64_t>>
using bit_vector_simd = bit_vector<allocator_t>;

}  // namespace libjst::utility
