// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the serialiser function for the jst.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/filesystem>

#include <cereal/archives/binary.hpp>

namespace jstmap
{

template <typename jst_t>
void serialise_jst(jst_t const & jst, cereal::BinaryOutputArchive & binary_archive)
{
    jst.save(binary_archive);
}

template <typename partitioned_jst_t>
void serialise_partitioned_jst(partitioned_jst_t const & partitioned_jst, cereal::BinaryOutputArchive & binary_archive)
{
    partitioned_jst.save(binary_archive);
}

template <typename jst_t, typename partitioned_jst_t>
void serialise(jst_t const & jst, partitioned_jst_t const & partitioned_jst,
               std::filesystem::path const & output_path)
{
    std::ofstream output_stream{output_path.c_str()};
    cereal::BinaryOutputArchive binary_archive{output_stream};
    serialise_jst(jst, binary_archive);
    serialise_partitioned_jst(partitioned_jst, binary_archive);
}

} // namespace jstmap
