// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <ranges>

#include <libio/file/formatted_file.hpp>
#include <libio/format/vcf/vcf_format.hpp>
#include <libio/format/vcf/vcf_record.hpp>

// TEST(vcf_file_test, vcf_header)
// {
//     // using file_t = libio::formatted_file<libio::vcf_record, libio::vcf_format>;
//     std::ifstream fstream{DATADIR"test_file.vcf"}; // opened fstream.
//     libio::tokenizer_streambuffer_adaptor buf{fstream.rdbuf()};
//     libio::vcf_token_header header_tkn{libio::vcf_token_meta{fstream}, libio::vcf_token_sample{fstream}};
//     // libio::pivot_tokenizer tknizer{buf, libio::pivot_matcher{"#CHROM"}};

//     // { // consume until done
//     //      libio::consume_tokenizer{std::move(tknizer)};
//     // }
//     // libio::line_tokenizer ltkn{buf};


//     // auto it = ltkn.begin();
//     // while (it != ltkn.end())
//     // {
//     //     ++it;
//     // }

//     // libio::vcf_token_meta meta_tkn{fstream};
//     libio::vcf_header header{};
//     libio::detokenize_to(header_tkn, header);

//     std::cout << "VER: " << header.version() << "\n";
//     std::cout << "INF: " << header.infos() << "\n";
//     std::cout << "SMP: " << header.sample_names() << "\n";

//     // std::cout << "version: " << header.version() << "\n";
//     // std::cout << "infos: " << header.infos() << "\n";
//     // // std::cout << "header: " << header.infos() << "\n";

// }

TEST(vcf_file_test, vcf_file)
{
    using file_t = libio::formatted_file<libio::vcf_record, libio::vcf_format>;
    file_t file{DATADIR"test_file.vcf", libio::vcf_format{}};
    auto const & fmt = file.format();
    std::cout << "Number of samples " << fmt.header().sample_names() << "\n";

    for (auto it = file.begin(); it != file.end(); ++it)
    {
        std::cout << "\n"
                  << "    CHROM: " << it->chrom() << "\n"
                  << "      POS: " << it->pos() << "\n"
                  << "       ID: " << it->id() << "\n"
                  << "      REF: " << it->ref() << "\n"
                  << "      ALT: " << it->alt() << "\n"
                  << "     QUAL: " << it->qual() << "\n"
                  << "   FILTER: " << it->filter() << "\n"
                  << "     INFO: " << it->info() << "\n"
                  << "   FORMAT: " << it->genotype_format() << "\n"
                  << "GENOTYPES: " << it->genotypes() << "\n";
    }
}

