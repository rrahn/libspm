// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides glue code to make it work with seqan q-gram index.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan/basic.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace seqan
{

template <>
struct ValueSize<seqan3::dna5>
{
    using Type = uint8_t;
    static constexpr Type VALUE = seqan3::alphabet_size<seqan3::dna5>;
};
} // namespace seqan

namespace seqan3
{
inline typename seqan::ValueSize<seqan3::dna5>::Type ordValue(seqan3::dna5 const & c) noexcept
{
    return seqan3::to_rank(c);
}
} // namespace seqan3
// ValueSize<TValue>::VALUE
// BitsPerValue<TValue>::VALUE?
// ordValue()
// plus cast?
