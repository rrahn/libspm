// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the function to parse a vcf file and construct a JST from it.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <atomic>

#include <seqan3/core/debug_stream/detail/to_string.hpp>

#include <jstmap/global/jstmap_type_alias.hpp>

namespace jstmap
{

// TODO:
//  * Make synchronised!
//  * Allow to print to file.
//!\brief An error handler to work with possible errors during file parsing.
class application_logger
{
private:

    bool _throw_on_error{true}; //!< Wether to throw or output a log message.
    verbosity_level _verbosity_level{verbosity_level::standard}; //!< Level of printed information in non-throwing mode.

public:

    constexpr application_logger() = default;
    constexpr application_logger(bool const throw_on_error, verbosity_level const level) noexcept :
        _throw_on_error{throw_on_error},
        _verbosity_level{level}
    {}

    //!\brief Logs the given message depending on the logger settings.
    template <typename ...message_args_t>
    constexpr void operator()(verbosity_level verbosity, logging_level log_level, message_args_t && ...message_args)
    {
        using namespace std::literals;

        bool const will_throw = log_level == logging_level::error && _throw_on_error;
        // Ignore this message if verbosity of message is not enabled by user.
        if (!will_throw && (verbosity == verbosity_level::quite ||
            static_cast<int>(verbosity) > static_cast<int>(_verbosity_level)))
            return;

        std::string_view msg_level{};

        switch (log_level)
        {
            case logging_level::info: msg_level = "[INFO] "sv; break;
            case logging_level::warning: msg_level = "[WARNING] "sv; break;
            case logging_level::error: msg_level = "[ERROR] "sv; break;
            default: throw std::runtime_error{"[ERROR] Unknown logging level!"};
        };

        auto message = seqan3::detail::to_string(msg_level, std::forward<message_args_t>(message_args)...);

        // If throw_on_error is enabled and the message indicates an error.
        if (will_throw)
            throw std::runtime_error{message};
        else
            std::cerr << message << "\n";
    }
};

/*!\brief Sets the new application wide logger.
 *
 * \param[in] handle The new logger handle to set.
 *
 * \returns The logger that was used before (never nullptr).
 *
 * \details
 *
 * Sets the new application logger. If the handle is a nullptr the logger handle is set
 * to the default application logger, which throws on error and uses a standard verbosity level.
 */
application_logger * set_application_logger(application_logger * handle) noexcept;

//!\brief Returns the application wide logger.
application_logger & get_application_logger() noexcept;

}  // namespace jstmap

