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

#include <exception>

#include <libcontrib/execute/concept_operation.hpp>
#include <libcontrib/execute/concept_sender.hpp>

namespace execute
{
    namespace _run
    {

        inline constexpr struct _fn
        {

            struct receiver
            {
                std::exception_ptr & _error;

                void set_value() const noexcept
                {}

                void set_done() const noexcept
                {}

                void set_error(std::exception_ptr error) noexcept
                {
                    _error = error;
                }

            };

            template <typename sender_t>
            void operator()(sender_t && sender) const {
                std::exception_ptr error{};
                auto run_command = execute::connect(sender, receiver{error});
                execute::start(run_command);

                if (error) {
                    std::rethrow_exception(error);
                }
            }
        } run;
    } // namespace _fn

    using _run::run;
} // namespace execute
