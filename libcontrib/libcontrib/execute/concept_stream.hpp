// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides stream concept.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libcontrib/std/tag_invoke.hpp>

namespace execute
{

    namespace _next {
        inline constexpr struct _cpo  {
            template <typename stream_t>
                requires std::tag_invocable<_cpo, stream_t>
            constexpr auto operator()(stream_t && stream) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, stream_t>)
                -> std::tag_invoke_result_t<_cpo, stream_t>
            {
                return std::tag_invoke(_cpo{}, (stream_t &&) stream);
            }

        private:

            template <typename stream_t>
                requires requires (stream_t stream) { { std::forward<stream_t &&>(stream).next() }; }
            constexpr friend auto tag_invoke(_cpo, stream_t && stream)
                noexcept(noexcept(std::declval<stream_t &&>().next()))
                -> decltype(std::declval<stream_t &&>().next())
            {
                return ((stream_t &&)stream).next();
            }
        } next;
    } // namespace _next
    using _next::next;

    template <typename stream_t>
    using next_t = std::tag_invoke_result_t<std::tag_t<execute::next>, stream_t>;

    namespace _cleanup {
        inline constexpr struct _cpo  {
            template <typename stream_t>
                requires std::tag_invocable<_cpo, stream_t>
            constexpr auto operator()(stream_t && stream) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo, stream_t>)
                -> std::tag_invoke_result_t<_cpo, stream_t>
            {
                return std::tag_invoke(_cpo{}, (stream_t &&) stream);
            }

        private:

            template <typename stream_t>
                requires requires (stream_t stream) { { std::forward<stream_t &&>(stream).cleanup() }; }
            constexpr friend auto tag_invoke(_cpo, stream_t && stream)
                noexcept(noexcept(std::declval<stream_t &&>().cleanup()))
                -> decltype(std::declval<stream_t &&>().cleanup())
            {
                return ((stream_t &&)stream).cleanup();
            }
        } cleanup;
    } // namespace _cleanup
    using _cleanup::cleanup;

    template <typename stream_t>
    using cleanup_t = std::tag_invoke_result_t<std::tag_t<execute::cleanup>, stream_t>;
}  // namespace execute
