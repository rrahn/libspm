// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>
#include <string>

#include <cereal/archives/binary.hpp>

#include <jstmap/search/load_jst.hpp>

namespace jstmap
{

std::pair<jst_t, std::unique_ptr<partitioned_jst_t>> load_jst(std::filesystem::path const & jst_input_file_path)
{
    using namespace std::literals;

    std::fstream jst_input_stream{jst_input_file_path};

    if (!jst_input_stream.good())
        throw std::runtime_error{"Couldn't open path for loading the jst! The path is ["s +
                                 jst_input_file_path.string() +
                                 "]"s};


    std::pair<jst_t, std::unique_ptr<partitioned_jst_t>> jst_handle{};
    cereal::BinaryInputArchive input_archive{jst_input_stream};
    jst_handle.first.load(input_archive);
    jst_handle.second = std::make_unique<partitioned_jst_t>(std::addressof(jst_handle.first));
    jst_handle.second->load(input_archive);

    return jst_handle;
}

} // namespace jstmap
