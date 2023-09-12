// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::utility::bit_vector_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <bit>
#include <compare>
#include <concepts>
#include <initializer_list>
#include <vector>

#include <cereal/types/base_class.hpp>

#include <seqan3/utility/detail/bits_of.hpp>

namespace libjst
{


/*!\brief A crtp-base class to provide generic access to the basic functions.
 *
 * \details
 *
 * This is an abstract base class which can only be instantiated thorugh its derived class.
 * The behaviour of the binary bit vector operations must be implemented by the derived class.
 */
template <typename derived_t, typename allocator_t>
class bit_vector_base :
    public std::vector<uint64_t, typename std::allocator_traits<allocator_t>::rebind_alloc<uint64_t>>
{
private:
    //!\brief The type of the underlying chunk of bits.
    using chunk_type = uint64_t;
    //!\brief The base type.
    using base_t = std::vector<chunk_type, typename std::allocator_traits<allocator_t>::rebind_alloc<chunk_type>>;

    //!\brief Befriend the derived class.
    friend derived_t;

    template <bool is_const>
    class bit_reference;

    template <bool is_const>
    class bit_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The iterator over the bits.
    using iterator = bit_iterator<false>;
    //!\brief The const iterator over the bits.
    using const_iterator = bit_iterator<true>;
    //!\brief The value type is `bool`.
    using value_type = std::iter_value_t<iterator>;
    //!\brief The reference type which is implemented as a proxy.
    using reference = std::iter_reference_t<iterator>;
    //!\brief The const reference type which is implemented as a proxy.
    using const_reference = std::iter_reference_t<const_iterator>;
    //!\brief The size_type.
    using size_type = size_t;
    //!\brief The difference type.
    using difference_type = std::iter_difference_t<iterator>;
    //!\brief The allocator type to use.
    using allocator_type = allocator_t;
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

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief The default constructor which optionally sets the allocator.
    bit_vector_base(allocator_t const & alloc = allocator_t{}) : base_t{alloc}
    {}

    bit_vector_base(bit_vector_base const &) = default; //!< Default.
    bit_vector_base(bit_vector_base &&) = default; //!< Default.
    bit_vector_base & operator=(bit_vector_base const &) = default; //!< Default.
    bit_vector_base & operator=(bit_vector_base &&) = default; //!< Default.
    ~bit_vector_base() = default; //!< Default.
    //!\}

public:
    /*!\name Member functions
     * \{
     */

    /*!\brief Assigns values to the container.
     *
     * \tparam iterator_t The type of the iterator; must model std::input_iterator.
     * \tparam sentinel_t The type of the sentinel; must model std::sentinel_for `iterator_t`.
     *
     * \param[in] first An first element to copy elements from.
     * \param[in] last The end of the range to copy elements from.
     *
     * \details
     *
     * Replaces the contents with copies of the range `[first, last)'. The behaviour is undefined if either argument
     * is an iterator to `*this`.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in distance between first and last.
     */
    template <std::input_iterator iterator_t, std::sentinel_for<iterator_t> sentinel_t>
    //!\cond
        requires std::assignable_from<bool &, std::iter_reference_t<iterator_t>>
    //!\endcond
    constexpr void assign(iterator_t first, sentinel_t last)
    {
        derived_t tmp{}; // To ensure strong exception guarantee.
        if constexpr (std::sized_sentinel_for<sentinel_t, iterator_t>)
            tmp.reserve(std::ranges::distance(first, last));

        std::ranges::copy(first, last, std::back_inserter(tmp));

        // ----- no exception after this.
        swap(tmp);
        set_new_size(std::ranges::distance(begin(), end()));
    }

    /*!\brief Assigns values to the container.
     *
     * \param[in] ilist The initialiser list with the elements to insert.
     *
     * \details
     *
     * Replaces the contents with the elements from the initializer list ilist.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in `ilist.size()`.
     */
    constexpr void assign(std::initializer_list<bool> const & ilist)
    {
        assign(std::ranges::begin(ilist), std::ranges::end(ilist));
    }

    /*!\brief Assigns values to the container.
     *
     * \param[in] count The new size of the container.
     * \param[in] bit The value to initialize elements of the container with.
     *
     * \details
     *
     * Replaces the contents with `count` copies of value `bit`.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in count.
     */
    constexpr void assign(size_type const count, bool const bit)
    {
        resize(count, bit);
        std::ranges::for_each(*as_base(), [value = fill_chunk(bit)] (chunk_type & chunk) { chunk = value;});
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Access specified element.
    constexpr reference operator[](difference_type const position) noexcept
    {
        assert(position >= 0);
        assert(static_cast<size_type>(position) < size());

        return *std::ranges::next(begin(), position);
    }

    //!\brief Access specified element.
    constexpr const_reference operator[](difference_type const position) const noexcept
    {
        assert(position >= 0);
        assert(static_cast<size_type>(position) < size());

        return *std::ranges::next(begin(), position);
    }

    /*!\brief Access the last element.
     *
     * \returns A reference to the last element in the container.
     *
     * \details
     *
     * Calling back on an empty container causes underfined behaviour.
     *
     * ### Exception
     *
     * Throws nothing.
     *
     * ### Complexity
     *
     * Constant.
     */
    constexpr reference back() noexcept
    {
        assert(!empty()); // Calling on empty container is undefined behaviour.

        return (*this)[size() - 1];
    }

    //!\overload
    constexpr const_reference back() const noexcept
    {
        assert(!empty()); // Calling on empty container is undefined behaviour.

        return (*this)[size() - 1];
    }

    //!\brief Checks if all bits are set to `true`.
    constexpr bool all() const noexcept
    {
        constexpr chunk_type mask = ~static_cast<chunk_type>(0);
        return std::ranges::all_of(*as_base(), [mask] (chunk_type const & chunk) { return chunk == mask; });
    }

    //!\brief Checks if any bit is set to `true`.
    constexpr bool any() const noexcept
    {
        constexpr chunk_type mask = static_cast<chunk_type>(0);
        return std::ranges::any_of(*as_base(), [mask] (chunk_type const & chunk) { return chunk | mask; });
    }

    //!\brief Checks if none of the bits is set to `true`.
    constexpr bool none() const noexcept
    {
        return !any();
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns the number of elements.
    constexpr size_type size() const noexcept
    {
        return _size;
    }

    //!\brief Checks wether the container is empty.
    constexpr bool empty() const noexcept
    {
        return _size == 0;
    }

    //!\brief Returns the capacity.
    constexpr size_type capacity() const noexcept
    {
        return base_t::capacity() * chunk_size;
    }

    /*!\brief Reserves storage.
     *
     * \param[in] new_capacity The new capacity of the bit vector.
     *
     * \details
     *
     * Increase the capacity of the vector to a value that's greater or equal to new_capacity. If new_capacity is
     * greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
     * reserve() does not change the size of the vector. If new_capacity is greater than capacity(), all iterators,
     * including the past-the-end iterator, and all references to the elements are invalidated. Otherwise, no
     * iterators or references are invalidated.
     *
     * ### Exceptions
     *
     * std::length_error if `new_capacity > max_size()` or any exception thrown by allocator_t::allocate().
     * If an exception is thrown this function has no effect (strong exception guarantee).
     */
    constexpr void reserve(size_type const new_capacity)
    {
        base_t::reserve(as_derived()->host_size_impl(new_capacity));
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Adds an element to the end.
     *
     * \param bit The bit to add to the end.
     *
     * \details
     *
     * Appends the given element value to the end of the container.
     *
     * If the new size() is greater than capacity() then all iterators and references
     * (including the past-the-end iterator) are invalidated.
     * Otherwise only the past-the-end iterator is invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown (which can be due to allocator_t::allocate(), this function has no effect
     * (strong exception guarantee).
     *
     * ### Complexity
     *
     * Amortised constant.
     */
    constexpr void push_back(bool bit)
    {
        size_t const new_size = size() + 1;
        resize(new_size);
        // ---- no exception after this point.
        set_new_size(new_size);
        back() = bit; // set the bit.
    }

    //!\brief Changes the number of elements stored, where additional copies of `bit` are appended.
    constexpr void resize(size_type const count, bool const bit = {})
    {
        base_t::resize(as_derived()->host_size_impl(count));

        size_t const old_size = size();
        set_new_size(count);
        if (size() > old_size) // If bit is true and we increase the size.
        {
            if (bit) std::ranges::fill(begin() + old_size, end(), bit);
        }
        else if (size() < old_size)
        {
            size_t const chunk_position = to_chunk_position(size());
            (*as_base())[chunk_position] &= (1 << to_local_chunk_position(size())) - 1;
            std::ranges::fill(std::ranges::next(as_base()->begin(), chunk_position + 1, as_base()->end()),
                              as_base()->end(),
                              0);
        }
    }

    //!\brief Performs binary AND between `this` and `rhs`.
    constexpr derived_t & operator&=(derived_t const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return as_derived()->binary_transform_impl(rhs, [] (auto const & left_chunk, auto const & right_chunk)
        {
            return left_chunk & right_chunk;
        });
    }

    //!\brief Performs binary OR between `this` and `rhs`.
    constexpr derived_t & operator|=(derived_t const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return as_derived()->binary_transform_impl(rhs, [] (auto const & left_chunk, auto const & right_chunk)
        {
            return left_chunk | right_chunk;
        });
    }

    //!\brief Performs binary XOR between `this` and `rhs`.
    constexpr derived_t & operator^=(derived_t const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return as_derived()->binary_transform_impl(rhs, [] (auto const & left_chunk, auto const & right_chunk)
        {
            return left_chunk ^ right_chunk;
        });
    }

    //!\brief Performs binary NOT.
    constexpr derived_t operator~() const noexcept
    {
        derived_t tmp{*as_derived()};
        return tmp.flip();
    }

    //!\brief Performs binary AND.
    constexpr friend derived_t operator&(derived_t lhs, derived_t const & rhs) noexcept
    {
        return lhs &= rhs;
    }

    //!\brief Performs binary OR.
    constexpr friend derived_t operator|(derived_t lhs, derived_t const & rhs) noexcept
    {
        return lhs |= rhs;
    }

    //!\brief Performs binary XOR.
    constexpr friend derived_t operator^(derived_t lhs, derived_t const & rhs) noexcept
    {
        return lhs ^= rhs;
    }

    //!\brief Computes the bitwise `a &= ~b` operator without an additional copy.
    constexpr derived_t & and_not(derived_t const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return as_derived()->binary_transform_impl(rhs, [] (auto const & left_chunk, auto const & right_chunk)
        {
            return left_chunk & ~right_chunk;
        });
    }

    //!\brief Flips all bits in-place.
    constexpr derived_t & flip() noexcept
    {
        std::ranges::for_each(*as_base(), [] (chunk_type & chunk) { chunk = ~chunk; });
        return *as_derived();
    }

    //!\brief Flips the bit at the given position.
    constexpr derived_t & flip(size_type position)
    {
        using namespace std::literals;

        if (position >= size())
            throw std::out_of_range{"The given posisiton "s +  std::to_string(position) +
                                    " is out of the range [0, "s + std::to_string(size()) + ")!"s};

        (*this)[position].flip();
        return *as_derived();
    }

    //!\brief Exchanges the contents of the container with those of others.
    constexpr void swap(derived_t & other) noexcept
    {
        base_t::swap(*other.as_base());
        std::swap(_size, other._size);
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

    //!\copydoc libjst::bit_vector_base::begin
    constexpr const_iterator begin() const noexcept
    {
        return const_iterator{base_t::data()};
    }

    //!\copydoc libjst::bit_vector_base::begin
    constexpr const_iterator cbegin() const noexcept
    {
        return begin();
    }

    //!\brief Returns an iterator to the end.
    constexpr iterator end() noexcept
    {
        return begin() + size();
    }

    //!\copydoc libjst::bit_vector_base::end
    constexpr const_iterator end() const noexcept
    {
        return begin() + size();
    }

    //!\copydoc libjst::bit_vector_base::end
    constexpr const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

    /*!\name Serialisation
     * \{
     */
    /*!\brief Saves this bit vector to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model typename.
     *
     * \param[in, out] archive The archive to serialise this object to.
     */
    template <typename output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(cereal::base_class<base_t>(this), _size);
    }

    /*!\brief Loads this this bit vector from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model typename.
     *
     * \param[in, out] archive The archive to serialise this object from.
     */
    template <typename input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(cereal::base_class<base_t>(this), _size);
    }
    //!\}

private:

    //!\brief Sets the new size.
    void set_new_size(size_type const new_size) noexcept
    {
        _size = new_size;
    }

    //!\brief Casts `this` to its base class.
    base_t const * as_base() const noexcept
    {
        return static_cast<base_t const *>(this);
    }

    //!\overload
    base_t * as_base() noexcept
    {
        return static_cast<base_t *>(this);
    }

    //!\brief Casts `this` to its derived class.
    derived_t const * as_derived() const noexcept
    {
        return static_cast<derived_t const *>(this);
    }

    //!\overload
    derived_t * as_derived() noexcept
    {
        return static_cast<derived_t *>(this);
    }

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
template <typename derived_t, typename allocator_t>
template <bool is_const>
class bit_vector_base<derived_t, allocator_t>::bit_reference
{
private:
    //!\brief Befriend the bit vector so it can instantiate this proxy with a particular position.
    template <typename, typename>
    friend class bit_vector_base;

    //!\brief The const or non-const chunk type to be represented.
    using maybe_const_chunk_type = std::conditional_t<is_const, chunk_type const, chunk_type>;

    maybe_const_chunk_type * _chunk{}; //!< The underlying chunk.
    chunk_type _chunk_mask{}; //!< The mask selecting the bit for this chunk.

    /*!\brief Constructs the refernce with represented bit position within the container.
     *
     * \param[in] chunk A pointer to the chunk that contains the represented bit.
     * \param[in] local_chunk_position The position of the bit within the chunk.
     */
    constexpr bit_reference(maybe_const_chunk_type * chunk, size_type const local_chunk_position) noexcept :
        _chunk{chunk},
        _chunk_mask{static_cast<chunk_type>(1) << to_local_chunk_position(local_chunk_position)}
    {}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    bit_reference() = delete; //!< Deleted.
    bit_reference(bit_reference const & other) = default;
    bit_reference(bit_reference && other) = default;
    bit_reference & operator=(bit_reference const & other) = default;
    bit_reference & operator=(bit_reference && other) = default;

    /*!\brief Assigns a bit to the referenced bit.
     *
     * \param[in] bit The bit to set.
     */
    constexpr bit_reference & operator=(bool const bit) noexcept
    {
        bit ? set() : clear();
        return *this;
    }

    //!\overload
    // Needed to model std::output_iterator<bit_iterator, bool>, which requires the assignment to an const && version
    // of the proxy.
    constexpr bit_reference const & operator=(bool const bit) const noexcept
    //!\cond
        requires (!is_const)
    //!\endcond
    {
        bit ? set() : clear();
        return *this;
    }
    //!\}

    //!\brief Converts this proxy to a bool.
    constexpr operator bool() const noexcept
    {
        return *_chunk & _chunk_mask;
    }

    //!\brief Flips the referenced bit.
    constexpr bit_reference & flip() noexcept
    {
        (*this) ? clear() : set();
        return *this;
    }

private:
    //!\brief Sets the bit at the specific position.
    constexpr void set() noexcept
    {
        *_chunk |= _chunk_mask;
    }

    //!\brief Clears the bit at the specific position.
    constexpr void clear() noexcept
    {
        *_chunk &= ~_chunk_mask;
    }
};

/*!\brief A random access iterator over the bit vector.
 * \implements std::random_access_iterator
 *
 * \tparam is_const A bool that indicates a const iterator if the value is `true`, or a non-const iterator otherwise.
 */
template <typename derived_t, typename allocator_t>
template <bool is_const>
class bit_vector_base<derived_t, allocator_t>::bit_iterator
{
private:
    //!\brief Befriend the bit_iterator types with different constness.
    template <bool>
    friend class bit_iterator;

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

    /*!\brief Copies from a non-const iterator.
     *
     * \param[in] other The other non-const iterator to copy from.
     */
    constexpr bit_iterator(bit_iterator<!is_const> const & other) noexcept requires (is_const) :
        _chunk{other._chunk},
        _chunk_position{other._chunk_position}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the currently pointer-to element.
    constexpr reference operator*() const noexcept
    {
        return reference{_chunk, _chunk_position};
    }

    //!\brief Returns the element `count` positions away from the current iterator position.
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

    //!\brief Returns a new iterator advanced by `count` many elements.
    constexpr bit_iterator operator+(difference_type const count) const noexcept
    {
        bit_iterator tmp{*this};
        return tmp += count;
    }

    //!\brief Returns a new iterator advanced by `count` many elements.
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

    //!\brief Advances the iterator by `count` many elements.
    constexpr bit_iterator & operator-=(difference_type const count) noexcept
    {
        return *this += -count;
    }

    //!\brief Returns a new iterator advances by `count` many elements.
    constexpr bit_iterator operator-(difference_type const count) const noexcept
    {
        bit_iterator tmp{*this};
        return tmp -= count;
    }

    //!\brief Returns the distance between `this` and the `rhs` iterator.
    template <bool is_const_other>
    constexpr difference_type operator-(bit_iterator<is_const_other> rhs) const noexcept
    {
        return ((_chunk - rhs._chunk) << division_mask) - // number of bits between chunks.
               to_local_chunk_position(rhs._chunk_position) + // minus the first bits in rhs.
               to_local_chunk_position(_chunk_position); // plus the first bits of the lhs
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares with another iterator.
    template <bool is_const_other>
    bool operator==(bit_iterator<is_const_other> const & rhs) const
    {
        return _chunk == rhs._chunk &&
              (to_local_chunk_position(_chunk_position) == to_local_chunk_position(rhs._chunk_position));
    }

    //!\brief Compares the two iterator by their chunk position and local chunk position.
    template <bool is_const_other>
    std::strong_ordering operator<=>(bit_iterator<is_const_other> const & rhs) const
    {
        if (std::strong_ordering order = _chunk <=> rhs._chunk; order == std::strong_ordering::equivalent)
            return to_local_chunk_position(_chunk_position) <=> to_local_chunk_position(rhs._chunk_position);
        else
            return order;
    }
    //!\}
};

}  // namespace libjst::utility
