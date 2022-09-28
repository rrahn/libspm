// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>
#include <ranges>

#include <seqan3/io/sequence_file/input.hpp>

#include <jstmap/search/load_queries.hpp>

namespace jstmap
{

struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = jst::contrib::dna5;
    using sequence_legal_alphabet = jst::contrib::dna5;
};

std::vector<raw_sequence_t> load_queries(std::filesystem::path const & query_input_file_path)
{
    seqan3::sequence_file_input<input_traits> query_input_file{query_input_file_path};

    std::vector<raw_sequence_t> queries{};

    std::ranges::for_each(query_input_file, [&] (auto const & sequence_record)
    {
        queries.emplace_back(seqan3::get<seqan3::field::seq>(sequence_record));
    });

    return queries;
}

} // namespace jstmap
