// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/execute/concept_receiver.hpp>

namespace execute
{
    namespace _ready_done
    {
        template <typename receiver_t>
        struct command
        {
            receiver_t _receiver;

            void start() noexcept
            {
                execute::set_done((receiver_t &&) _receiver);
            }
        };

        struct sender
        {
            template <typename receiver_t>
            command<receiver_t> connect(receiver_t &&receiver) const &
            {
                return command<receiver_t>{(receiver_t &&) receiver};
            }
        };
    } // namespace _ready_done

    using ready_done_sender = _ready_done::sender;
} // namespace execute
