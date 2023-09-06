// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides an application logger to record different messages.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <jstmap/global/application_logger.hpp>

namespace jstmap
{

namespace
{

//!\brief The default logger of this application.
constinit application_logger default_application_logger{};
//!\brief The application handle which handles the logger in a sequentially synchronized fashion.
constinit std::atomic<application_logger *> application_logger_handle{&default_application_logger};

} // namespace

//!\cond
// See documentation of declaration in <jstmap/index/application_logger.hpp>
application_logger * set_application_logger(application_logger * logger) noexcept
{
    if (logger == nullptr)
        logger = &default_application_logger;

    return application_logger_handle.exchange(logger);
}

application_logger & get_application_logger() noexcept
{
    return *application_logger_handle.load();
}
//!\endcond

}  // namespace jstmap
