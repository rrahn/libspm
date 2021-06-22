// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap::test
{

using seqan3::operator""_dna5;

inline const raw_sequence_t reference = "aacctt"_dna5;

inline const std::vector<raw_sequence_t> sequences
{
    "aaaaaa"_dna5,
    "cccccc"_dna5,
    "tttttt"_dna5
};

}  // namespace jstmap::test
