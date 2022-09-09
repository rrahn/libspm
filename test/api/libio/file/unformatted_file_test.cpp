// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <libio/format/fasta/fasta_format.hpp>
#include <libio/file/unformatted_file.hpp>

TEST(unformatted_file_test, basic)
{
    using namespace std::literals;

    libio::unformatted_file file{DATADIR"in.fasta", libio::fasta_format{}};
    auto record = file.read_record();
    ASSERT_EQ(record.seq(), "ACGTTTGATTCGCG"s);
    ASSERT_EQ(record.id(), "seq1"s);
}

TEST(unformatted_file_test, valid_extensions)
{
    ASSERT_EQ(libio::valid_extensions(libio::fasta_format{}),
              (std::vector<std::string>{".fa", ".fasta", ".fn"}));
}

TEST(unformatted_file_test, formatted)
{
    using namespace std::literals;

    libio::unformatted_file file{DATADIR"in.fasta", libio::fasta_format{}};

    // Idea: add transformation layer!
    // - this reinterprets the records into the user specified types
    // - the transformation is handled by a CPO call which the user can overload for the given container type
    // - we can offer seqan3 specific transformations and also extend it with accelerated versions
    // - this can be easily extended.
    // - we may also wrap then a given record in order to reduce the stuff that needs to be converted
    // - finally we use this abstraction to wrap around a range-based adaptor interface.
        // - steps: load buffer from the stream block wise
        // - in each block find the delimiter pattern
        // - store the byte offset for each delimter
        // - can be reused in multiple places

    auto record = file.read_record();
    ASSERT_EQ(record.seq(), "ACGTTTGATTCGCG"s);
    ASSERT_EQ(record.id(), "seq1"s);
}
