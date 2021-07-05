// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functions to write the results.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <vector>

#include <seqan/sequence.h>

#include <seqan3/io/sam_file/output.hpp>

#include <libjst/context_position.hpp>

#include <jstmap/search/search_queries.hpp>

namespace jstmap
{

using sam_file_t = seqan3::sam_file_output<seqan3::fields<seqan3::field::ref_offset,
                                                          seqan3::field::seq,
                                                          seqan3::field::alignment>>;

void write_results(sam_file_t &,
                   std::vector<search_match> const &,
                   std::vector<raw_sequence_t> const &);

}  // namespace jstmap
