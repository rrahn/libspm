// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides serialiser object for the jst.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/std/tag_invoke.hpp>
#include <libcontrib/type_traits.hpp>

#include <libjst/container/concept_serialiser.hpp>

namespace libjst
{

    template <typename archive_t, typename target_t>
    class delegate_serialiser_impl
    {
        // we might need to get the value to be pointed to?
        archive_t _archive;
        std::reference_wrapper<target_t> _target;

    public:
        template <typename _archive_t, typename _target_t>
        delegate_serialiser_impl(_archive_t &&archive, _target_t &target) :
            _archive{(_archive_t &&) archive},
            _target{target}
        {
        }

        // every one is a normal archiver -> we pass through to the wrapped archiver!
        template <typename... args_t>
        void operator()(args_t &&...args)
        {
            _archive((args_t &&) args...);
        }

        template <typename external_target_t>
        requires std::same_as<external_target_t, std::remove_cvref_t<target_t>>
        constexpr friend auto tag_invoke(std::tag_t<libjst::load_extern>,
                                         delegate_serialiser_impl &me,
                                         [[maybe_unused]] external_target_t const &external_target)
        {
            assert(std::addressof(me._target.get()) == std::addressof(external_target)); // reading the same object
            libjst::load(me._target.get(), me._archive);                                 // delegate the load process
        }

        template <typename external_target_t>
        requires std::same_as<external_target_t, std::remove_cvref_t<target_t>>
        constexpr friend auto tag_invoke(std::tag_t<libjst::save_extern>,
                                         delegate_serialiser_impl &me,
                                         external_target_t const &external_target)
        {
            assert(std::addressof(me._target.get()) == std::addressof(external_target)); // reading the same object
            libjst::save(me._target.get(), me._archive);
        }
    };

    template <typename archive_t, typename target_t>
    delegate_serialiser_impl(archive_t &&, target_t &) -> delegate_serialiser_impl<archive_t, target_t>;

    namespace _delegate_serialiser
    {
        inline constexpr struct _fn
        {
            template <typename serialiser_t, typename target_t>
            constexpr auto operator()(serialiser_t &&serialiser, target_t &target) const noexcept
            {
                using original_t = jst::contrib::maybe_unwrap_t<target_t>;
                return libjst::delegate_serialiser_impl{(serialiser_t &&) serialiser, static_cast<original_t &>(target)};
            }

            template <typename target_t>
            constexpr auto operator()(target_t &target) const noexcept
                -> jst::contrib::closure_result_t<_fn, std::reference_wrapper<target_t>>
            {
                return jst::contrib::make_closure(_fn{}, std::ref(target));
            }
        } delegate_serialiser;
    } // namespace _delegate_serialiser
    using _delegate_serialiser::delegate_serialiser;

} // namespace libjst
