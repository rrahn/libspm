// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides search base algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <exception>

#include <libcontrib/std/tag_invoke.hpp>

#include <libjst/traversal/concept.hpp>

namespace libjst
{

    // ----------------------------------------------------------------------------
    // Algorithm CPOs
    // ----------------------------------------------------------------------------

    // what default can we use here with the algorithm?
    namespace _default_search
    {
        template <typename delegate_t>
        struct runnable
        {
            delegate_t _delegate;
            std::exception_ptr & _error;

        private:
            template <typename this_t, typename result_t>
                requires std::same_as<std::remove_cvref_t<this_t>, runnable>
            friend void tag_invoke(std::tag_t<libjst::set_next>, this_t && me, result_t && result) noexcept
            { // publish result
                std::invoke(me._delegate, (result_t &&)result);
            }

            template <typename this_t, typename error_t>
                requires std::same_as<std::remove_cvref_t<this_t>, runnable>
            friend void tag_invoke(std::tag_t<libjst::set_error>, this_t && me, error_t const & error) noexcept
            { // store error.
                me._error = std::make_exception_ptr(error);
            }

            template <typename cpo_t, typename this_t>
                requires std::same_as<std::remove_cvref_t<this_t>, runnable>
            friend void tag_invoke(cpo_t, this_t &&) noexcept
            { // noop
            }
        };

        inline constexpr struct _cpo_t {
            template <typename haystack_t, typename searcher_t, typename callback_t>
                requires sender_to<std::invoke_result_t<searcher_t, haystack_t>, runnable<std::decay_t<callback_t>>>
            constexpr void operator()(haystack_t &&haystack, searcher_t &&searcher, callback_t && callback) const
            {
                std::exception_ptr error;
                auto op = libjst::connect(searcher(haystack), runnable<std::decay_t<callback_t>>{(callback_t &&)callback, error});
                libjst::start(op);

                if (error)
                    std::rethrow_exception(error);
            }
        } default_search;
    } // namespace _default_search
    using _default_search::default_search;

    namespace _search {
        inline constexpr struct _cpo_t {
            template <typename haystack_t, typename searcher_t, typename callback_t>
                requires std::tag_invocable<_cpo_t, haystack_t, searcher_t, callback_t>
            constexpr auto operator()(haystack_t &&haystack, searcher_t &&searcher, callback_t && callback) const
                noexcept(std::is_nothrow_tag_invocable_v<_cpo_t, haystack_t, searcher_t, callback_t>)
                -> std::tag_invoke_result_t<_cpo_t, haystack_t, searcher_t, callback_t>
            {
                return std::tag_invoke(_cpo_t{}, (haystack_t &&)haystack, (searcher_t &&)searcher, (callback_t &&)callback);
            }

            // defaults: try to invoke the algorithm directly, but we might get a sender back?
            template <typename haystack_t, typename searcher_t, typename callback_t>
                requires (!std::tag_invocable<_cpo_t, haystack_t, searcher_t, callback_t>)
            constexpr void operator()(haystack_t &&haystack, searcher_t &&searcher, callback_t &&callback) const
            {
                // now we can search for defaults!
                if constexpr (std::invocable<std::tag_t<default_search>, haystack_t, searcher_t, callback_t>) {
                    default_search((haystack_t &&)haystack, (searcher_t &&)searcher, (callback_t &&)callback);
                }
            }
        } search_base;
    }
    using _search::search_base;

}  // namespace libjst
