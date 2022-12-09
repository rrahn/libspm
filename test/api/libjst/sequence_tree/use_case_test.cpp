// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include "rcs_store_mock.hpp"

#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/sequence_node_tree.hpp>
#include <libjst/sequence_tree/covered_node_tree.hpp>

TEST(use_case_test, state_oblivious) {

    jst::test::mock_store<std::string> rcs_store{};

    // basic tree interface:
    libjst::volatile_tree base_tree{rcs_store};
    libjst::sequence_node_tree<decltype(base_tree)> labeled_tree{base_tree};
    libjst::covered_node_tree<decltype(labeled_tree)> covered_tree{labeled_tree};


    // steps to construct the jst
        // 1. make rcs_store a rooted store
        // 2. set the coverage information

    // move-operation and
    // data
    root(jst), sink(jst) // -> means wrappable
    node::label() // -> wrappable

    // volatile tree factories
    jst::tree::fragment_sequence_tree{rcs_store}; // each tree node represents a single sequence.
    jst::tree::root_path_sequence_tree{rcs_store}; // label sequence represents entire path from root to current node.

    // returns a node
    // the node itself can be the context?

    jst::tree::covered_sequence_variation_tree{rcs_store, initial_coverage};

    // communication within trees via context
    // how does tree access label?

    // tree adapter
    tree::transform(tree, fn_operation) {
        -> general purpose adapter.
        -> does not change the traversal
    }

    tree::reverse(tree) {
        -> new tree that traverses from sink to root
    }

    tree::trim(tree, window_size) {
        -> new tree that
    }

    tree::observe(tree, listener) {
        -> new tree that
    }

    tree::prune_uncovered(tree) {
        -> new tree that only handles covered trees
    }

    tree::chunk(tree, split_size) {
        -> returns a forest :) -> collection of trees.
    }

    // tree sinks
    // traversal independent
    tree::traverse_alt_first()
    // ...
    tree::traverse_flat()
}

