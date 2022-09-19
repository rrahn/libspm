// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <seqan3/std/ranges>
#include <sstream>

#include <libio/file/formatted_stream.hpp>
#include <libio/format/fasta/fasta_format.hpp>
#include <libio/format/fastq/fastq_format.hpp>
#include <libio/format/sequence/sequence_format.hpp>
#include <libio/format/sequence/sequence_record.hpp>

// Now we just want to tokenize it!

inline constexpr std::string_view fa_input = R"fa(>SEQ_ID 1
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
AGACTGAGCTACGAGCTAGCGACT
>SEQ_ID 2
GGTTAAGGTTCCCCAAGGTTAC
)fa";

inline constexpr std::string_view fq_input = R"fq(@SEQ_ID 1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID 2
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
)fq";

TEST(formatted_stream_test, fasta_format)
{
    std::stringstream sstream{fa_input.data()};
    libio::sequence_format seq_fmt{libio::fasta_format{}, libio::fastq_format{}};
    libio::select_format(seq_fmt, std::filesystem::path{"tmp.fa"});
    libio::formatted_stream fmt_stream{seq_fmt, std::move(sstream)};
    libio::sequence_record seq_record1{};
    libio::sequence_record seq_record2{};
    fmt_stream >> seq_record1;
    ASSERT_FALSE(fmt_stream.eof());
    fmt_stream >> seq_record2; // so we need to parse to the end of the record or the end of the stream?
    ASSERT_TRUE(fmt_stream.eof());

    std::cout << "First:\nID: " << seq_record1.id() << "\nSEQ: " << seq_record1.seq() << "\n\n";
    std::cout << "Second:\nID: " << seq_record2.id() << "\nSEQ: " << seq_record2.seq() << "\n\n";
    ASSERT_TRUE(fmt_stream.eof());
}

