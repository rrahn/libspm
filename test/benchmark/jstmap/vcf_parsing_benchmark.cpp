// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <vector>

#include <seqan/vcf_io.h>

#include <seqan3/test/performance/units.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/range/container/aligned_allocator.hpp>

#include <libio/file/formatted_file.hpp>
#include <libio/format/vcf/vcf_format.hpp>
#include <libio/format/vcf/vcf_record.hpp>

#include "benchmark_utility.hpp"

template <typename... args_t>
static void seqan2_vcf(benchmark::State &state, args_t &&...args)
{
    std::filesystem::path vcf_file_path{get<0>(std::tuple{args...})};

    seqan::VcfFileIn vcf_file{vcf_file_path.c_str()};
    seqan::VcfHeader vcf_header{};
    seqan::readHeader(vcf_header, vcf_file);

    seqan::VcfRecord _record{};

    std::vector<size_t> records_per_chr{};
    records_per_chr.resize(1000, 0);
    for (auto _ : state)
    {
        while (!seqan::atEnd(vcf_file))
        {
            seqan::readRecord(_record, vcf_file);
            ++records_per_chr[_record.rID];
        }
    }

    size_t vcf_file_size = file_size(vcf_file_path);
    size_t record_count{};
    std::ranges::for_each(records_per_chr, [&](size_t hits)
                          { record_count += hits; });

    state.counters["file_size"] = vcf_file_size;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(vcf_file_size);
    state.counters["#records"] = record_count;
}

class my_vcf_record : public libio::vcf_record
{
private:
    using base_t = libio::vcf_record;

public:
    my_vcf_record() = default;

private:
    template <typename buffer_t>
    constexpr friend auto tag_invoke(libio::tag_t<libio::set_field>,
                                     my_vcf_record &me,
                                     libio::field_code_type<libio::vcf_field::chrom> const &fc,
                                     buffer_t &&ibuffer) noexcept
    {
        return libio::set_field(static_cast<base_t &>(me), fc, (buffer_t &&) ibuffer);
    }

    template <libio::vcf_field field_tag, typename buffer_t>
        requires(field_tag != libio::vcf_field::chrom)
    constexpr friend auto tag_invoke(libio::tag_t<libio::set_field>,
                                     my_vcf_record &,
                                     libio::field_code_type<field_tag> const &,
                                     buffer_t &&) noexcept
    {
    }
};

template <typename... args_t>
static void libio_vcf(benchmark::State &state, args_t &&...args)
{
    std::filesystem::path vcf_file_path{get<0>(std::tuple{args...})};

    libio::vcf_format fmt{};
    libio::formatted_file<my_vcf_record, libio::vcf_format> vcf_file{vcf_file_path};

    std::vector<size_t> records_per_chr{};
    records_per_chr.resize(1000, 0);
    size_t copied_bytes{};
    for (auto _ : state)
    {
        // try {
        for (auto const &record : vcf_file)
        {
            ++records_per_chr[record.chrom()];
            copied_bytes += record.bytes();
        }
        // } catch (...) {
        //     size_t record_count{};
        //     std::ranges::for_each(records_per_chr, [&] (size_t hits) { record_count += hits; });
        //     std::cerr << "record_count = " << record_count << "\n";
        //     std::rethrow_exception(std::current_exception());
        // }
    }

    size_t vcf_file_size = file_size(vcf_file_path);
    size_t record_count{};
    std::ranges::for_each(records_per_chr, [&](size_t hits)
                          { record_count += hits; });

    state.counters["file size"] = vcf_file_size;
    state.counters["record count"] = record_count;
    state.counters["copied bytes"] = copied_bytes;
    state.counters["b_stream/s"] = seqan3::test::bytes_per_second(vcf_file_size);
    state.counters["#records/s"] = records_per_second(record_count);
    state.counters["reads/s"] =  seqan3::test::bytes_per_second(vcf_file_size + copied_bytes);
}

template <typename... args_t>
static void libio_vcf_stream(benchmark::State &state, args_t &&...args)
{
    std::filesystem::path vcf_file_path{get<0>(std::tuple{args...})};
    libio::vcf_format fmt{};
    libio::formatted_stream vcf_stream{fmt, std::ifstream{vcf_file_path.c_str()}};

    std::vector<char, seqan3::aligned_allocator<char, 64>> buffer{};
    buffer.resize(1024*1024);
    vcf_stream.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

    if (!vcf_stream.good())
        throw std::runtime_error{"no good stream."};

    using pos_type = typename std::remove_reference_t<decltype(vcf_stream)>::pos_type;

    std::vector<pos_type> token_positions{};
    token_positions.reserve(10000);

    for (auto _ : state)
    {
        while (!vcf_stream.eof())
        {
            { // load only tokens
                // std::cout << "before...   ";
                auto token = vcf_stream.get(); // no detokenization
                token_positions.push_back(token.position());
                // std::cout << "   ...after\n";
            }
        }
    }

    size_t vcf_file_size = file_size(vcf_file_path);

    state.counters["file size"] = vcf_file_size;
    state.counters["token cont"] = token_positions.size();
    state.counters["b_stream/s"] = seqan3::test::bytes_per_second(vcf_file_size);
}

// BENCHMARK_CAPTURE(seqan2_vcf, vcf, DATADIR "1KGP.chr22_20k.vcf");
BENCHMARK_CAPTURE(libio_vcf_stream, vcf, DATADIR"1KGP.chr22.vcf");
// BENCHMARK_CAPTURE(libio_vcf, vcf, DATADIR"1KGP.chr22.vcf");

BENCHMARK_MAIN();
