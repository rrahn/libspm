// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the jstmap wide application logger.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <atomic>
#include <iostream>

#include <seqan3/core/debug_stream/detail/to_string.hpp>

namespace jstmap
{

//!\brief An enum to select the verbosity level.
enum class verbosity_level : uint8_t
{
    quite, //!< No logging output is emitted.
    standard, //!< Logs regular information with no extra information on the output.
    verbose //!< Extra verbose logging output for debugging purposes.
};

//!\brief An enum to select the logging level.
enum class logging_level : uint8_t
{
    info, //!< An informative message during the execution.
    warning, //!< A warning message for non-severe issues during the execution.
    error, //!< An error message for severe issues during the execution.
    debug //!< A debug information only printed when verbose logging is enabled.
};

// TODO:
//  * Make synchronised!
//  * Allow to print to file.
//!\brief An error handler to work with possible errors during file parsing.
class application_logger
{
private:

    bool _throw_on_error{false}; //!< Wether to throw or output a log message.
    verbosity_level _verbosity_level{verbosity_level::standard}; //!< Level of printed information in non-throwing mode.

public:

    constexpr application_logger() = default;
    constexpr application_logger(bool const throw_on_error, verbosity_level const level) noexcept :
        _throw_on_error{throw_on_error},
        _verbosity_level{level}
    {}

    constexpr void set_verbosity(verbosity_level new_level) noexcept {
        _verbosity_level = new_level;
    }

    //!\brief Logs the given message depending on the logger settings.
    template <typename ...message_args_t>
    constexpr void operator()(verbosity_level verbosity, logging_level log_level, message_args_t && ...message_args)
    {
        using namespace std::literals;

        if ((_verbosity_level == verbosity_level::quite) ||
            (log_level == logging_level::debug && _verbosity_level != verbosity_level::verbose))
            return;

        std::string_view msg_level{};

        switch (log_level)
        {
            case logging_level::info: msg_level = "[INFO] "sv; break;
            case logging_level::warning: msg_level = "[WARNING] "sv; break;
            case logging_level::error: msg_level = "[ERROR] "sv; break;
            case logging_level::debug: msg_level = "[DEBUG] "sv; break;
            default: throw std::runtime_error{"[ERROR] Unknown logging level!"};
        };

        auto message = seqan3::detail::to_string(msg_level, std::forward<message_args_t>(message_args)...);
        // If throw_on_error is enabled and the message indicates an error.
        std::cerr << message << "\n";
    }

    template <typename ...message_args_t>
    constexpr void log(logging_level log_level, message_args_t && ...message_args) {
        this->operator()(_verbosity_level, log_level, (message_args_t &&) message_args...);
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

template <typename ...message_args_t>
constexpr void log(logging_level log_level, message_args_t && ...message_args) {
    get_application_logger().log(log_level, (message_args_t &&) message_args...);
}

template <typename ...message_args_t>
constexpr void log_info(message_args_t && ...message_args) {
    get_application_logger().log(logging_level::info, (message_args_t &&) message_args...);
}

template <typename ...message_args_t>
constexpr void log_debug(message_args_t && ...message_args) {
    get_application_logger().log(logging_level::debug, (message_args_t &&) message_args...);
}

template <typename ...message_args_t>
constexpr void log_warn(message_args_t && ...message_args) {
    get_application_logger().log(logging_level::warning, (message_args_t &&) message_args...);
}

template <typename ...message_args_t>
constexpr void log_err(message_args_t && ...message_args) {
    get_application_logger().log(logging_level::error, (message_args_t &&) message_args...);
}

}  // namespace jstmap

