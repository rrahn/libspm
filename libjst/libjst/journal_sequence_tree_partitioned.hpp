// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree partitioned.
 * \author Tom Lukas Lankenau <tom.lankenau AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/core/detail/strong_type.hpp>

#include <libjst/journal_sequence_tree_context_enumerator.hpp>
#include <libjst/journaled_sequence_tree.hpp>

namespace libjst
{

//!\brief A strong type to pass a context size object.
struct context_size : public seqan3::detail::strong_type<uint32_t, context_size>
{
    //!\brief The base type.
    using base_t = seqan3::detail::strong_type<uint32_t, context_size>;
    using base_t::base_t;
};

//!\brief A strong type to pass a bin index object.
struct bin_index : public seqan3::detail::strong_type<uint32_t, bin_index>
{
    //!\brief The base type.
    using base_t = seqan3::detail::strong_type<uint32_t, bin_index>;
    using base_t::base_t;
};

/*!\brief A partitioned journal sequence tree.
 *
 * \tparam jst_t The type of the journal sequence tree to wrap.
 *
 * \details
 *
 * This wrapper handles a collection of libjst::detail::journal_sequence_tree_traverser_model with non-overlapping
 * intervals over the given jst. Later the jst can be traversed in bins by constructing the respective agent for a
 * particular bin. The partitioned jst is a wrapper class and cannot be default initialised.
 */
template <typename jst_t>
class journal_sequence_tree_partitioned
{
public:
    //!\brief The type of the traverser model to manage.
    using traverser_model_t = detail::journal_sequence_tree_traverser_model<jst_t>;
    //!\brief The type of the context enumerator to use.
    using context_enumerator_type = typename jst_t::context_enumerator_type;

private:
    std::vector<traverser_model_t> _bins{}; //!\< The container stroing the model for each bin.
    jst_t const * _jst{}; //!\< A shared pointer to the underlying jst.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_partitioned() = delete; //!< Delete.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned const &) = default; //!< Default.
    constexpr journal_sequence_tree_partitioned(journal_sequence_tree_partitioned &&) = default; //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_partitioned & operator=(journal_sequence_tree_partitioned &&)
        = default; //!< Default.
    ~journal_sequence_tree_partitioned() = default; //!< Default.

    /*!\brief Constructs the partitioned jst from a jst pointer and a bin count.
     *
     * \param[in] jst A pointer to the jst to wrap.
     * \param[in] bin_count The number of bins. Defaults to 1 and must be `> 0`.
     */
    journal_sequence_tree_partitioned (jst_t const * jst, size_t bin_count = 1) : _jst{jst}
    {
       assert(_jst != nullptr);
       assert(bin_count > 0);

        std::ptrdiff_t bin_size = (_jst->reference().size() + bin_count - 1) / bin_count;

        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(bin_count); i++)
            _bins.push_back(traverser_model_t{_jst, i * bin_size, (i + 1) * bin_size});
    }
    //!\}

    //!\brief Returns the bin count.
    size_t bin_count() const noexcept
    {
        return _bins.size();
    }

    //!\brief Returns a new context enumerator from the given bin and the context size.
    context_enumerator_type context_enumerator(libjst::context_size const context_size,
                                               libjst::bin_index const bin_index) const
    {
        if (bin_index.get() >= bin_count())
            throw std::out_of_range{"The bin index: " + std::to_string(bin_index.get()) + " is out of range!"};

        return context_enumerator_type{_bins[bin_index.get()], context_size.get()};
    }
};

} // namespace libjst
