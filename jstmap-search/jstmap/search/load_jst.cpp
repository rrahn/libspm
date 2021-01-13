// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/just_map/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>

#include <cereal/archives/binary.hpp>

#include <jstmap/search/load_jst.hpp>

namespace jstmap
{

jst_t load_jst(std::filesystem::path const & jst_input_file_path)
{
    std::fstream jst_input_stream{jst_input_file_path};

    jst_t jst{};

    {
        cereal::BinaryInputArchive input_archive{jst_input_stream};
        jst.load(input_archive);
    }

    return jst;
}

} // namespace jstmap
