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
#include <libio/format/fasta/fasta_format.hpp>
#include <libio/format/fastq/fastq_format.hpp>
#include <libio/format/sequence/sequence_format.hpp>
#include <libio/format/sequence/sequence_record.hpp>

TEST(formatted_file_test, fasta_format)
{
    libio::sequence_format seq_fmt{libio::fasta_format{}, libio::fastq_format{}};
    using file_t = libio::formatted_file<libio::sequence_record, decltype(seq_fmt)>;
    file_t file{DATADIR"in.fasta", std::move(seq_fmt)};
    for (auto it = file.begin(); it != file.end(); ++it)
    {
        std::cout << "ID: " << it->id() << "\nSEQ: " << it->seq() << "\n\n";
    }
}

