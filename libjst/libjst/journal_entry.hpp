// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <set>
#include <type_traits>
#include <ranges>

namespace libjst::detail
{

// insertion sequence is different then reference sequence
template <std::ranges::view segment_t>
//!\cond
    requires std::ranges::sized_range<segment_t>
//!\endcond
class journal_entry
{
private:

    mutable segment_t _segment{};
    mutable uint32_t _segment_begin_position{};

public:

    using size_type = uint32_t;
    using segment_type = segment_t;

    journal_entry() = default;
    journal_entry(journal_entry const &) = default;
    journal_entry(journal_entry &&) = default;
    journal_entry & operator=(journal_entry const &) = default;
    journal_entry & operator=(journal_entry &&) = default;
    ~journal_entry() = default;

    journal_entry(size_type segment_begin_position, segment_type segment) :
        _segment{std::move(segment)},
        _segment_begin_position{segment_begin_position}
    {}

    // size_type & segment_begin_position() noexcept
    // {
    //     return _segment_begin_position;
    // }

    size_type & segment_begin_position() const noexcept
    {
        return _segment_begin_position;
    }

    size_type segment_end_position() const noexcept
    {
        return segment_begin_position() + segment_size();
    }

    // segment_type & segment() noexcept
    // {
    //     return _segment;
    // }

    segment_type & segment() const noexcept //noexcept(std::is_nothrow_copy_constructible_v<segment_type>)
    {
        return _segment;
    }

    size_type segment_size() const noexcept
    {
        return std::ranges::size(_segment);
    }

    constexpr bool operator==(journal_entry const &) const = default;

    constexpr std::weak_ordering operator<=>(journal_entry const & rhs) const noexcept
    {
        return segment_end_position() <=> rhs.segment_end_position();
    }

    constexpr std::weak_ordering operator<=>(size_type const & position) const noexcept
    {
        return segment_end_position() <=> position;
    }
};

template <typename segment_t>
journal_entry(uint32_t, segment_t) -> journal_entry<segment_t>;

template <typename char_t, typename char_traits_t, typename segment_t>
std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream,
                                                       journal_entry<segment_t> const & je)
{
    stream << "[" << je.segment_begin_position() << ", " << je.segment_end_position() << ") : ";
    stream << "<";
    for (auto c : je.segment())
        stream << c;
    stream << ">";

    return stream;
}
}  // namespace libjst::detail
