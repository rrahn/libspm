// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector_adaptor.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <compare>
#include <concepts>
#include <vector>

#include <seqan3/std/bit>
#include <seqan3/utility/detail/bits_of.hpp>

namespace libjst
{

template <typename t>
concept allocator = requires (std::allocator_traits<t> alloc_traits)
{
    typename decltype(alloc_traits)::value_type;
};

/*!\brief An allocator aware bit vector.
 *
 * \details
 *
 * Implements an [allocator aware](https://en.cppreference.com/w/cpp/named_req/AllocatorAwareContainer) bit vector
 * on the basis of a std::vector with `uint64_t` as value type. The bit vector can be dynamically resized and
 * provides additional interfaces to apply efficient bit-operations on it.
 *
 * The reference type is a special proxy that provides access to a single bit. Note that it is not a real reference
 * but can be converter to a bool or assigned from a bool.
 */
template <typename allocator_t = std::allocator<bool>>
class bit_vector : public std::vector<uint64_t, typename std::allocator_traits<allocator_t>::rebind_alloc<uint64_t>>
{
private:
    //!\brief The type of the underlying chunk of bits.
    using chunk_type = uint64_t;
    //!\brief The base type.
    using base_t = std::vector<chunk_type, typename std::allocator_traits<allocator_t>::rebind_alloc<chunk_type>>;

    template <bool is_const>
    class bit_reference;

    template <bool is_const>
    class bit_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type is `bool`.
    using value_type = bool;
    //!\brief The reference type which is implemented as a proxy.
    using reference = bit_reference<false>;
    //!\brief The const reference type which is implemented as a proxy.
    using const_reference = bit_reference<true>;
    //!\brief The size_type.
    using size_type = size_t;
    //!\brief The iterator over the bits.
    using iterator = bit_iterator<false>;
    //!\brief The const iterator over the bits.
    using const_iterator = bit_iterator<true>;
    //!\}

private:
    // ----------------------------------------------------------------------------
    // static constexpr member variables
    // ----------------------------------------------------------------------------

    //!\brief The number of bits represented in one chunk, e.g. 64.
    static constexpr size_type chunk_size = seqan3::detail::bits_of<chunk_type>;
    //!\brief The mask used for the modulo operations using the bitwise and operator, e.g. & 63.
    static constexpr size_type modulo_mask = chunk_size - 1;
    //!\brief The mask used for the division operations using bitwise shift operator, e.g. >> 6.
    static constexpr size_type division_mask = std::countr_zero(chunk_size);

    // ----------------------------------------------------------------------------
    // member variables
    // ----------------------------------------------------------------------------

    size_type _size{}; //!< The number of elements.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    bit_vector() = default; //!< Default.

    /*!\brief Constructs the bit vector with `count` copies of elements with value `bit`.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] bit The bit to set during initialisation.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector(size_type const count, bool const bit, allocator_t const & alloc = allocator_t{}) :
        base_t{chunks_needed(count), fill_chunk(bit), alloc},
        _size{count}
    {}

    /*!\brief Constructs the container with `count` default-inserted instances of `bool`. No copies are made.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] alloc The allocator to use [optional].
     */
    constexpr bit_vector(size_type const count, allocator_t const & alloc = allocator_t{}) :
        bit_vector{count, bool{}, alloc}
    {}
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns the number of elements.
    constexpr size_type size() const noexcept
    {
        return _size;
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the beginning.
    constexpr iterator begin() noexcept
    {
        return iterator{base_t::data()};
    }

    //!\copydoc libjst::bit_vector::begin
    constexpr const_iterator begin() const noexcept
    {
        return const_iterator{base_t::data()};
    }

    //!\copydoc libjst::bit_vector::begin
    constexpr const_iterator cbegin() const noexcept
    {
        return begin();
    }

    //!\brief Returns an iterator to the end.
    constexpr iterator end() noexcept
    {
        return begin() + size();
    }

    //!\copydoc libjst::bit_vector::end
    constexpr const_iterator end() const noexcept
    {
        return begin() + size();
    }

    //!\copydoc libjst::bit_vector::end
    constexpr const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

private:
    //!\brief Returns how many chunks are needed to store `count` many elements.
    constexpr size_type chunks_needed(size_type const count) const noexcept
    {
        return (count + 63) >> 6; // ceil(count/64)
    }

    //!\brief Returns a new chunk filled with the given bit.
    constexpr chunk_type fill_chunk(bool const bit) const noexcept
    {
        return (bit) ? ~static_cast<chunk_type>(0) : 0;
    }

    //!\brief Converts the position to the local position within the chunk.
    static constexpr size_type to_local_chunk_position(size_type const position) noexcept
    {
        return position & modulo_mask; // e.g. position % 64
    }

    //!\brief Converts the position to the chunk position.
    static constexpr size_type to_chunk_position(size_type const position) noexcept
    {
        return position >> division_mask; // e.g. position / 64
    }
};

/*!\brief The bit proxy returned as reference.
 *
 * \tparam is_const A bool that indicates a const proxy if the value is `true`, or a non-const proxy otherwise.
 *
 * \details
 *
 * This proxy is returned as a proxy for the reference type since a single bit cannot be addressed directly.
 * This proxy allows all operations that can be used on a single bit and is also implicitly convertible to a bool.
 * It cannot be default constructed and can only be instantiated with a particular bit from the bit vector and its
 * associated classes.
 */
template <typename allocator_t>
template <bool is_const>
class bit_vector<allocator_t>::bit_reference
{
private:
    //!\brief Befriend the bit vector so it can instantiate this proxy with a particular position.
    template <typename>
    friend class bit_vector;

    //!\brief The const or non-const chunk type to be represented.
    using maybe_const_chunk_type = std::conditional_t<is_const, chunk_type const, chunk_type>;

    maybe_const_chunk_type * _chunk{}; //!< The underlying chunk.
    size_type _chunk_position{}; //!< The bit position within the chunk.

    /*!\brief Constructs the refernce with represented bit position within the container.
     *
     * \param[in] chunk A pointer to the chunk that contains the represented bit.
     * \param[in] local_chunk_position The position of the bit within the chunk.
     */
    constexpr bit_reference(maybe_const_chunk_type * chunk, size_type const local_chunk_position) noexcept :
        _chunk{chunk},
        _chunk_position{to_local_chunk_position(local_chunk_position)}
    {}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    bit_reference() = delete; //!< Deleted.
    //!\}

    //!\brief Converts this proxy to a bool.
    constexpr operator bool() const noexcept
    {
        return *_chunk & (1 << _chunk_position);
    }
};

/*!\brief A random access iterator over the bit vector.
 * \implements std::random_access_iterator
 *
 * \tparam is_const A bool that indicates a const iterator if the value is `true`, or a non-const iterator otherwise.
 */
template <typename allocator_t>
template <bool is_const>
class bit_vector<allocator_t>::bit_iterator
{
private:
    //!\brief The type of the chunk.
    using maybe_const_chunk_type = std::conditional_t<is_const, chunk_type const, chunk_type>;

    maybe_const_chunk_type * _chunk{}; //!< The underlying chunk.
    size_type _chunk_position{}; //!< The bit position within the chunk.

public:
    /*!\name Associated types
     * \{
     */
    using value_type = bool; //!< The value type.
    using reference = bit_reference<is_const>; //!< The proxy type used as reference.
    using pointer = void; //!\< The pointer type is void.
    using difference_type = std::ptrdiff_t; //!< The difference type.
    using iterator_category = std::random_access_iterator_tag; //!< The iterator category.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bit_iterator() = default; //!< Default.

    /*!\brief Constructs the iterator set to the begin of the given chunk.
     *
     * \param[in] chunk A pointer to the chunk that contains the represented bit.
     */
    explicit constexpr bit_iterator(maybe_const_chunk_type * chunk) noexcept : _chunk{chunk}, _chunk_position{0}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the currently pointer-to bit.
    constexpr reference operator*() const noexcept
    {
        return reference{_chunk, _chunk_position};
    }

    //!\brief Returns the bit at `count` position right from the current iterator position.
    constexpr reference operator[](difference_type const count) const noexcept
    {
        return *((*this) + count);
    }
    //!\}

    /*!\name Arithmetic operator
     * \{
     */
    //!\brief Increments the iterator by one.
    constexpr bit_iterator & operator++() noexcept
    {
        _chunk += !static_cast<bool>(to_local_chunk_position(++_chunk_position));
        return *this;
    }

    //!\brief Increments the iterator by one and returns the iterator before the increment.
    constexpr bit_iterator operator++(int) noexcept
    {
        bit_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    //!\brief Advances the iterator by `count` many elements.
    constexpr bit_iterator & operator+=(difference_type const count) noexcept
    {
        //           chunk:|    0   |    1   |    2   |    3   |    4   |    5   |
        //                 |--------|--------|--------|--------|-x------|--------|
        //  chunk_position:|01234567|01234567|01234567|01234567|01234567|01234567|
        // global position:|01234567|89012345|67890123|45678901|23456789|01234567|
        //                 |0       |  1     |    2   |      3 |        |4       |
        if (count < 0)
        {
            size_type updated_count = modulo_mask - to_local_chunk_position(_chunk_position) - count;
            _chunk_position = modulo_mask - to_local_chunk_position(updated_count);
            _chunk -= to_chunk_position(updated_count); //(to_chunk_position(-count) + (old_chunk_position < _chunk_position));
        }
        else
        {
            _chunk += to_chunk_position(to_local_chunk_position(_chunk_position) + count);
            _chunk_position = to_local_chunk_position(_chunk_position + count);
        }

        return *this;
    }

    //!\brief Returns a new iterator incremented by `count` many elements.
    constexpr bit_iterator operator+(difference_type const count) const noexcept
    {
        bit_iterator tmp{*this};
        return tmp += count;
    }

    //!\brief Returns a new iterator incremented by `count` many elements.
    friend constexpr bit_iterator operator+(difference_type const count, bit_iterator rhs) noexcept
    {
        return rhs + count;
    }

    //!\brief Decrements the iterator by one.
    constexpr bit_iterator & operator--() noexcept
    {
        _chunk -= !static_cast<bool>(to_local_chunk_position(--_chunk_position));
        return *this;
    }

    //!\brief Decrements the iterator by one and returns the iterator before the decrement.
    constexpr bit_iterator operator--(int) noexcept
    {
        bit_iterator tmp{*this};
        --(*this);
        return tmp;
    }

    //!\brief Decrements the iterator by `count` many elements.
    constexpr bit_iterator & operator-=(difference_type const count) noexcept
    {
        return *this += -count;
    }

    //!\brief Returns a new iterator decremented by `count` many elements.
    constexpr bit_iterator operator-(difference_type const count) const noexcept
    {
        bit_iterator tmp{*this};
        return tmp -= count;
    }

    //!\brief Returns the distance between this and the `rhs` iterator.
    constexpr difference_type operator-(bit_iterator rhs) const noexcept
    {
        return ((_chunk - rhs._chunk) << division_mask) - // number of bits between chunks.
               to_local_chunk_position(rhs._chunk_position) + // minus the first bits in rhs.
               to_local_chunk_position(_chunk_position); // plus the first bits of the lhs
    }
    //!\}

    /*!\name Comparison operator
     * \{
     */
    //!\brief Compares with another iterator.
    bool operator==(bit_iterator const & rhs) const
    {
        return _chunk == rhs._chunk &&
              (to_local_chunk_position(_chunk_position) == to_local_chunk_position(rhs._chunk_position));
    }

    //!\brief Compares the two iterator by their chunk position and local chunk position.
    std::strong_ordering operator<=>(bit_iterator const & rhs) const
    {
        if (std::strong_ordering order = _chunk <=> rhs._chunk; order == std::strong_ordering::equivalent)
            return to_local_chunk_position(_chunk_position) <=> to_local_chunk_position(rhs._chunk_position);
        else
            return order;
    }
    //!\}
};

}  // namespace libjst::utility
