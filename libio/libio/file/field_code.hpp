// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides field_code implementation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    // class field_code; // forward declaration

    // class field_category
    // {
    // public:
    //     constexpr field_category() noexcept = default;
    //     field_category(const & field_category) = delete; // non-copyable/movable/assignable
    //     virtual ~field_category() = default;
    //     virtual std::string_view name() const noexcept = 0;

    //     // these are virtual overloads!
    //     // can we make them overload differently?
    //     // Basically, this is a CPO now?
    //     virtual bool equivalent(int32_t fc, field_condition const & fcond)  const noexcept
    //     {
    //         return false;
    //     }

    //     virtual bool equivalent(field_code const & fc, int fcond)  const noexcept
    //     {
    //         return false;
    //     };

    // };// abstract base class?
    // this to customize by the user for the specific format/record
    // something that returns a field_code?
    template <auto field_tag>
    class field_code_type
    {

        // TODO how can we set the field category?
        // field_category const * _category;

    public:

        using value_type = decltype(field_tag);

        static constexpr value_type value = field_tag;

        constexpr field_code_type() noexcept
        {
            // find overload for the field tag?
            // requirements for field_tag?
            // select_field_category(field_tag)
        }

        // template <typename field_category_t>
        //     // requires is_field_category_v<field_category_t>
        // constexpr implicit field_code(int32_t fid, field_category_t const & fcat) :
        //     _field_id{fid},
        //     _category{std::addressof(fcat)}
        // {}

        // template <typename field_code_t>
        //     // requires is_field_code_v<field_id_t> // not a real concept?
        // constexpr implicit field_code(field_code_t fc) : field_code{libio::make_field_code(fc)}
        // {}

        constexpr value_type operator()() const noexcept
        {
            return value;
        }

        // constexpr field_category const & category() const noexcept
        // {
        //     return *_category;
        // }

        constexpr operator bool() const noexcept
        {
            return value != 0; // if field_id != 0 then it is set.
        }
    };

    template <auto field_tag>
    inline constexpr field_code_type<field_tag> field_code{};
    // overload is_error_code_enum? can we do this more liniently?
    // need to plug a unique global object
    //

    // usage scenario:
    // - register enum value as field code
    // - define equivalence between different field codes? // error_condition & error_code
    // - define domain/category of field_code
    // - need to be templatized? we use it to decide if we can find an alternative field selction?
    // - wether the user wants to provide its own customised access to the records, is not our business.

} // namespace libio
