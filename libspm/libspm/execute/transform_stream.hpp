// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides make stream factory.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libspm/execute/ready_done.hpp>
#include <libspm/execute/then.hpp>

#include <libspm/copyable_box.hpp>
#include <libspm/execute/concept_stream.hpp>

namespace execute
{

    namespace _transform_stream
    {

        template <typename parent_stream_t, typename fn_t>
        class stream {

            using fn_box_t = jst::contrib::copyable_box<std::remove_reference_t<fn_t>>;

            parent_stream_t _parent_stream;
            fn_box_t _fn;

        public:

            explicit stream(parent_stream_t parent_stream, fn_t fn) :
                _parent_stream{(parent_stream_t&&) parent_stream},
                _fn{(fn_t&&)fn}
            {}

            auto next() noexcept {
                return execute::next(_parent_stream) | execute::then(std::ref(_fn));
            }

            auto cleanup() noexcept {
                return execute::cleanup(_parent_stream);
            }
        };

        inline struct closure
        {
            template <typename parent_stream_t, typename fn_t>
            auto operator()(parent_stream_t&& parent_stream, fn_t&& fn) const
                noexcept(std::is_nothrow_constructible_v<stream<parent_stream_t, fn_t>>)
                -> stream<parent_stream_t, fn_t>
            {
                static_assert(std::is_rvalue_reference_v<fn_t &&>);
                return stream<parent_stream_t, fn_t>{(parent_stream_t &&) parent_stream, (fn_t &&) fn};
            }

            template <typename fn_t>
            auto operator()(fn_t && fn) const
                noexcept(noexcept(jst::contrib::make_closure(std::declval<closure>(), (fn_t&&)fn)))
                -> jst::contrib::closure_result_t<closure, fn_t>
            {
                static_assert(std::is_rvalue_reference_v<fn_t &&>);
                return jst::contrib::make_closure(closure{}, (fn_t &&)fn);
            }

        } transform_stream;

    } // namespace _transform_stream

    using _transform_stream::transform_stream;
} // namespace execute
