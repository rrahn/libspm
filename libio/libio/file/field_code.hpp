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

#include <concepts>
#include <type_traits>

#include <libio/utility/tag_invoke.hpp>

namespace libio
{
    template <typename field_tag_t>
    struct field_code_category;

    template <typename field_tag_t>
        // requres is_format_field_tag_v<>
    using field_code_category_t = typename field_code_category<field_tag_t>::type;

    // this to customize by the user for the specific format/record
    // something that returns a field_code?
    template <auto field_tag>

    class field_code_type
    {
        // TODO how can we set the field category?
        // field_category const * _category;
        // we could set the category and then translate to something differently?
        using category_t = field_code_category_t<decltype(field_tag)>; // more on a static level.
        // so we can set the category based on the field set.
        // we want to compare?
        // what do we want to do?
        // and what can we do?
        category_t _field_category{}; // do we need it to be constexpr?

    public:

        using value_type = decltype(field_tag);

        static constexpr value_type value = field_tag;

        constexpr field_code_type() noexcept
        {
        }

        // it is already set?
        // constexpr implicit field_code(int32_t fid, category_t const & fcat) :
        //     _field_id{fid},
        //     _category{std::addressof(fcat)}
        // {}

        // template <typename field_tag_t>
        //     // requires is_field_code_v<field_id_t> // not a real concept?
        // constexpr implicit field_code(field_tag_t field_tag) : field_code{libio::make_field_code(fc)}
        // {}

        constexpr category_t const & category() const noexcept
        {
            return _field_category;
        }

        constexpr value_type operator()() const noexcept
        {
            return value;
        }

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

    namespace _equivalent
    {
        inline constexpr struct _cpo
        {
            // can we not have a default that switches the arguments?
            // independent of category the sam
            template <typename category_t, auto field_tag, typename field_tag_t>
                requires (std::same_as<decltype(field_tag), field_tag_t>)
            constexpr friend bool tag_invoke(_cpo,
                                             category_t const &,
                                             field_code_type<field_tag> const & fc,
                                             field_tag_t const & ft) noexcept
            {
                return fc() == ft;
            }

            template <typename category_t, auto field_tag, typename field_tag_t>
                requires tag_invocable<_cpo, category_t const &, field_code_type<field_tag> const &, field_tag_t const &>
            constexpr bool operator()(category_t const & category,
                                      field_code_type<field_tag> const & fc,
                                      field_tag_t const & ft) const noexcept
            {
                return libio::tag_invoke(_cpo{}, category, fc, ft);
            }

            template <typename category_t, typename field_tag_t, auto field_tag>
                requires tag_invocable<_cpo, category_t const &, field_code_type<field_tag> const &, field_tag_t const &>
            constexpr bool operator()(category_t const & category,
                                      field_tag_t const & ft,
                                      field_code_type<field_tag> const & fc) const noexcept
            {
                return (*this)(category, fc, ft);
            }
        } equivalent{};
    }
    using _equivalent::equivalent;

    class default_field_code_category  // base implementation
    {
    public:

        constexpr default_field_code_category() = default;

    private:

        template <auto field_tag, typename field_tag_t>
            requires (!std::same_as<decltype(field_tag), field_tag_t>)
        constexpr friend bool tag_invoke(tag_t<libio::equivalent>,
                                         default_field_code_category const &,
                                         field_code_type<field_tag> const &,
                                         field_tag_t const &) noexcept
        {
            return false;
        }
    };

        // requres is_format_field_tag_v<>
    template <typename field_tag_t>
    struct field_code_category : public std::type_identity<default_field_code_category>
    {
    };

} // namespace libio
