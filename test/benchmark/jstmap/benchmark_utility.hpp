// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <filesystem>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#include <jstmap/create/vcf_parser.hpp>
#include <jstmap/global/application_logger.hpp>
#include <jstmap/global/load_jst.hpp>

using sequence_t = seqan3::dna5_vector;

auto create_jst_from_vcf(std::filesystem::path reference_file, std::filesystem::path vcf_file)
{
    jstmap::application_logger logger{false, jstmap::verbosity_level::quite};
    jstmap::set_application_logger(&logger);
    return std::move(jstmap::construct_jst_from_vcf(reference_file, vcf_file).front());
}

template <typename sequences_t, typename algorithm_t>
auto naive_traversal(sequences_t && sequences, algorithm_t && algorithm)
{
    std::ranges::for_each(sequences, [&] (auto const & haystack)
    {
        auto haystack_range = std::ranges::subrange{std::ranges::begin(haystack), std::ranges::end(haystack)};
        while (!std::ranges::empty(haystack_range))
            haystack_range = algorithm(haystack_range);
    });
}

sequence_t generate_query(size_t const query_size)
{
    using alphabet_t = std::ranges::range_value_t<sequence_t>;
    return seqan3::test::generate_sequence<alphabet_t>(query_size);
}
