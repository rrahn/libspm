// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>
#include <string>

#include <cereal/archives/binary.hpp>

#include <jstmap/global/load_jst.hpp>

namespace jstmap
{

rcs_store_t load_jst(std::filesystem::path const & rcs_store_path)
{
    using namespace std::literals;

    std::fstream rcs_store_stream{rcs_store_path};
    if (!rcs_store_stream.good())
        throw std::runtime_error{"Couldn't open path for loading the jst! The path is ["s +
                                 rcs_store_path.string() +
                                 "]"s};

    rcs_store_t rcs_store{};
    {
        cereal::BinaryInputArchive input_archive{rcs_store_stream};
        rcs_store.load(input_archive);
    }
    return rcs_store;
}

} // namespace jstmap
