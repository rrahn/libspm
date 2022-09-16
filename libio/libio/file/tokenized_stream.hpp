// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of the tokenized strea.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iostream>

namespace libio
{
template <typename token_t>
class tokenized_stream : protected std::iostream {
private:

    using base_t = std::iostream;

    // We can now manage a special input/output iterator that directly works on the buffer set by the iostream

public:

    tokenized_stream() = default;
    tokenized_stream(std::basic_streambuf<CharT,Traits>* stream_buffer) : base_t{stream_buffer}
    {
    }

    token_t get() {
        return token_t{base_t::rdbuf()}; // we initialise the token with the rdbuf()
        // Implement the tokenized parsing via the iterator over the virtual stream.
        // here we create a token from the stream.
    }

    void put(token_t && token) {
        // we can use the tokenisation function to return a token?
        // here we write the token to the stream.
    }

    template <typename record_t>
    void put(record_t && record) {
        token_t tmp{base_t::rdbuf()};
        libio::tokenize(tmp, std::forward<record_t>(record)); // now tokenize the record which dumps the record into the buffer?
    }

private:
};
}  // namespace libio
