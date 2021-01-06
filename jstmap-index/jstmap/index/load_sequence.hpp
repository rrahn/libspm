// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/filesystem>
#include <vector>

#include <jstmap/index/global_types.hpp>

namespace jstmap
{

std::vector<raw_sequence_t> load_sequences(std::filesystem::path const &);

} // namespace jstmap
