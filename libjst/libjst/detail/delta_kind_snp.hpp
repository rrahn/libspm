// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::delta_kind_snp.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <span>
#include <vector>

#include <cereal/types/base_class.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/range/sequence.hpp>

#include <libjst/detail/delta_kind_base.hpp>

namespace libjst::detail
{

/*!\brief A delta event representing a substitution.
 *
 * \tparam alphabet_t The alphabet type used to store the inserted sequence; must model seqan3::semialphabet.
 */
template <seqan3::semialphabet alphabet_t>
class delta_kind_snp : public delta_kind_base<seqan3::alphabet_rank_t<alphabet_t>>
{
private:
    using base_t = delta_kind_base<seqan3::alphabet_rank_t<alphabet_t>>; //!< The base type.

    // requires alphabet to be constexpr constructible?!
    static const inline std::array<alphabet_t, seqan3::alphabet_size<alphabet_t>> values
    {
        [] () constexpr
        {
            using rank_t = std::remove_cvref_t<decltype(seqan3::alphabet_size<alphabet_t>)>;
            std::array<alphabet_t, seqan3::alphabet_size<alphabet_t>> tmp{};
            for (rank_t idx = 0; idx < tmp.size(); ++idx)
                seqan3::assign_rank_to(idx, tmp[idx]);

            return tmp;
        } ()
    };

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr delta_kind_snp() = default; //!< Default.
    constexpr delta_kind_snp(delta_kind_snp const &) = default; //!< Default.
    constexpr delta_kind_snp(delta_kind_snp &&) = default; //!< Default.
    constexpr delta_kind_snp & operator=(delta_kind_snp const &) = default; //!< Default.
    constexpr delta_kind_snp & operator=(delta_kind_snp &&) = default; //!< Default.
    ~delta_kind_snp() = default; //!< Default.

    /*!\brief Initialises a substitution with the given substitution sequence.
     *
     * \tparam sequence_t The sequence type; must model seqan3::sequence.
     *
     * \param[in] sequence The substituted sequence.
     */
    template <seqan3::sequence sequence_t>
    explicit constexpr delta_kind_snp(sequence_t && sequence) : base_t{seqan3::to_rank(*std::ranges::begin(sequence))}
    {
        assert(std::ranges::distance(sequence) == 1);
    }
    //!\}

    //!\brief Returns the stored value.
    std::span<alphabet_t const, 1> value() const noexcept
    {
        using span_t = std::span<alphabet_t const, 1>;
        return span_t{values.begin() + static_cast<base_t const &>(*this).value(), 1};
    }


    //!\brief Compare against other substitutions for equality.
    bool operator==(delta_kind_snp const &) const = default;

    /*!\name Serialisation
     * \{
     */
    //!\copydoc libjst::detail::delta_kind_base::save
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(cereal::base_class<base_t>(this));
    }

    //!\copydoc libjst::detail::delta_kind_base::load
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(cereal::base_class<base_t>(this));
    }
    //!\}
};

/*!\name Type deduction guide
 * \{
 */
/*!\brief Deduces the alphabet type from the given sequence.
 * \relates libjst::detail::delta_kind_snp
 */
template <seqan3::sequence sequence_t>
delta_kind_snp(sequence_t) -> delta_kind_snp<std::ranges::range_value_t<sequence_t>>;
//!\}
}  // namespace libjst::detail
