// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief The method to save the index.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <jstmap/index/save_index.hpp>

namespace jstmap
{

void save_index(seqan3::interleaved_bloom_filter<> & ibf, index_options const & options)
{
    std::ofstream ostr{options.output_file};
    cereal::BinaryOutputArchive oarch{ostr};
    oarch(options.bin_size);
    ibf.serialize(oarch);
}

} // namespace jstmap
