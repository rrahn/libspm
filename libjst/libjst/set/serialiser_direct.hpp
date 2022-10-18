// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides direct serialiser class.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <libcontrib/closure_adaptor.hpp>
#include <libcontrib/std/tag_invoke.hpp>
#include <libcontrib/type_traits.hpp>

#include <libjst/set/concept_serialiser.hpp>

namespace libjst
{
    template <typename archive_t, typename value_t>
    class direct_serialiser_impl
    {
        // we might need to get the value to be pointed to?
        archive_t _archive;
        std::reference_wrapper<value_t> _value;

    public:
        template <typename _archive_t, typename _value_t>
        direct_serialiser_impl(_archive_t &&archive, _value_t &value) : _archive((_archive_t &&) archive), _value(value)
        {
        }

        template <typename... args_t>
        void operator()(args_t &&...args)
        {
            _archive((args_t &&) args...);
        }

    private:
        // don't need that here because it is a source.
        constexpr friend auto tag_invoke(std::tag_t<libjst::load_extern>,
                                         direct_serialiser_impl &me,
                                         [[maybe_unused]] value_t const &value)
        {
            assert(std::addressof(me._value.get()) == std::addressof(value)); // must point to the same memory location.
            me._archive(me._value.get());
        }

        constexpr friend auto tag_invoke(std::tag_t<libjst::save_extern>,
                                         direct_serialiser_impl &me,
                                         [[maybe_unused]] value_t const &value)
        {
            assert(std::addressof(me._value.get()) == std::addressof(value)); // must point to the same memory location.
            me._archive(me._value.get());
        }
    };

    template <typename archive_t, typename value_t>
    direct_serialiser_impl(archive_t &&, value_t &) -> direct_serialiser_impl<archive_t, value_t>;

    namespace _direct_serialiser
    {
        inline constexpr struct _fn
        {
            template <typename serialiser_t, typename target_t>
            constexpr auto operator()(serialiser_t &&serialiser, target_t &target) const noexcept
            {
                using original_t = jst::contrib::maybe_unwrap_t<target_t>;
                return libjst::direct_serialiser_impl{(serialiser_t &&) serialiser, static_cast<original_t &>(target)};
            }

            template <typename target_t>
            constexpr auto operator()(target_t &target) const noexcept
                -> jst::contrib::closure_result_t<_fn, std::reference_wrapper<target_t>>
            {
                return jst::contrib::make_closure(_fn{}, std::ref(target));
            }
        } direct_serialiser;
    } // namespace _direct_serialiser

    using _direct_serialiser::direct_serialiser;

}  // namespace libjst
