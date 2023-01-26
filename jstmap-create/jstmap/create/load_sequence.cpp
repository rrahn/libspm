// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>

#include <jstmap/create/load_sequence.hpp>

namespace jstmap
{

sequence_collection_t load_sequences(std::filesystem::path const & sequence_file)
{
    sequence_collection_t sequences{};

    seqan3::sequence_file_input<sequence_input_traits> fin{sequence_file.c_str()};
    std::ranges::for_each(fin, [&] (auto const & sequence_record)
    {
        sequences.push_back(sequence_record.sequence());
    });

    return sequences;
}

} // namespace jstmap
