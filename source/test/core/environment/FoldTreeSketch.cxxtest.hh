// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/environment/FoldTreeSketch.hh>

//Other headers
#include <core/kinematics/FoldTree.hh>

#include <core/types.hh>

#include <test/core/init_util.hh>

//C++ headers
#include <iostream>


class FoldTreeSketch : public CxxTest::TestSuite {
public:


  // Shared data elements go here.

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  // --------------- Test Cases --------------- //
  void test_simple_res_sketch(){
    core::environment::FoldTreeSketch fts;
    TS_ASSERT_THROWS_NOTHING( fts = core::environment::FoldTreeSketch( 1 ) );

    //Verify size.
    TS_ASSERT_EQUALS( fts.nres(), core::Size( 1 ) );

    //Invalid insertion calls
    TS_ASSERT_THROWS( fts.insert_cut( 0 ), core::environment::EXCN_FTSketchGraph );
    TS_ASSERT_THROWS( fts.insert_cut( 1 ), core::environment::EXCN_FTSketchGraph );
    TS_ASSERT_THROWS( fts.insert_jump( 1, 2 ), core::environment::EXCN_FTSketchGraph );

    //Invalid property queries
    TS_ASSERT( !fts.has_cut( 1 ) );
    TS_ASSERT_THROWS( fts.has_cut( 0 ), core::environment::EXCN_FTSketchGraph );
    TS_ASSERT_THROWS( fts.has_jump( 1, 2 ), core::environment::EXCN_FTSketchGraph );

    //Check rendering
    core::kinematics::FoldTreeOP ft_out;
    TS_ASSERT_THROWS_NOTHING( ft_out = fts.render() );
    TS_ASSERT( ft_out->is_simple_tree() );
    core::kinematics::FoldTree simple_ft;
    simple_ft.simple_tree( 1 );
    TS_ASSERT_EQUALS( *ft_out, simple_ft );

    //Ensure no cycles are detected
    utility::vector1< core::Size > cycle;
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 1 ) );
    TS_ASSERT_EQUALS( cycle.size(), core::Size( 0 ) );

    //Now with a bigger sketch
    TS_ASSERT_THROWS_NOTHING( fts = core::environment::FoldTreeSketch( 10 ) );
    TS_ASSERT_EQUALS( fts.nres(), core::Size( 10 ) );

    //Invalid insertions
    TS_ASSERT_THROWS( fts.insert_cut( 0 ), core::environment::EXCN_FTSketchGraph );
    TS_ASSERT_THROWS( fts.insert_cut( 10 ), core::environment::EXCN_FTSketchGraph );
    TS_ASSERT_THROWS( fts.insert_jump( 5, 11 ), core::environment::EXCN_FTSketchGraph );

    //Check graph integrity
    TS_ASSERT( !fts.has_cut( 1 ) );
    TS_ASSERT( !fts.has_cut( 9 ) );
    TS_ASSERT( !fts.has_jump( 1, 10 ) );
    TS_ASSERT( !fts.has_jump( 4, 7 ) );

    //Ensure no false cycles are detected.
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 1 ) );
    TS_ASSERT_EQUALS( cycle.size(), core::Size( 0 ) );

    //Check rendering
    TS_ASSERT_THROWS_NOTHING( fts.render( *ft_out ) );
    TS_ASSERT( ft_out->is_simple_tree() );
    simple_ft.simple_tree( 10 );

    TS_ASSERT_EQUALS( *ft_out, simple_ft );
  }

  void test_one_jump_res_sketch(){
    core::environment::FoldTreeSketch fts;
    TS_ASSERT_THROWS_NOTHING( fts = core::environment::FoldTreeSketch( 10 ) );

    //Test successful jump insertion
    TS_ASSERT_THROWS_NOTHING( fts.insert_jump( 4, 7 ) );
    TS_ASSERT( !fts.has_jump( 4, 6 ) );
    TS_ASSERT( fts.has_jump( 4, 7 ) );
    TS_ASSERT( !fts.has_jump( 4, 8 ) );

    //Verify cannot build cyclic graph
    TS_ASSERT_THROWS( fts.render(), core::environment::EXCN_FTSketchGraph );

    //Verify cycles are detected correctly.
    utility::vector1< core::Size > cycle;
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 1 ) );

    std::set< core::Size > cycle_set;
    std::copy( cycle.begin(), cycle.end(), std::inserter( cycle_set, cycle_set.end() ) );

    TS_ASSERT_DIFFERS( cycle_set.find( 4 ), cycle_set.end() );
    TS_ASSERT_DIFFERS( cycle_set.find( 5 ), cycle_set.end() );
    TS_ASSERT_DIFFERS( cycle_set.find( 6 ), cycle_set.end() );
    TS_ASSERT_DIFFERS( cycle_set.find( 7 ), cycle_set.end() );

    //Verify successful cut addition
    TS_ASSERT_THROWS_NOTHING( fts.insert_cut( 5 ) );
    TS_ASSERT( !fts.has_cut( 4 ) );
    TS_ASSERT( fts.has_cut( 5 ) );
    TS_ASSERT( !fts.has_cut( 6 ) );

    //Verify cycle is not detected.
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 7 ) );
    TS_ASSERT_EQUALS( cycle.size(), core::Size( 0 ) );

    //Verify successfully build graph
    core::kinematics::FoldTreeOP ft_out;
    TS_ASSERT_THROWS_NOTHING( ft_out = fts.render() );
    TS_ASSERT_EQUALS( fts.nres(), ft_out->nres() );
    TS_ASSERT( ft_out->is_cutpoint( 5 ) );
    TS_ASSERT_EQUALS( ft_out->num_cutpoint(), 1 );
    TS_ASSERT_EQUALS( ft_out->jump_nr( 4, 7 ), core::Size( 1 ) );
    TS_ASSERT_EQUALS( ft_out->num_jump(), core::Size( 1 ) );
  }

  void test_cycle_removal() {
    core::environment::FoldTreeSketch fts;
    TS_ASSERT_THROWS_NOTHING( fts = core::environment::FoldTreeSketch( 20 ) );

    // set up two-jump, no-cut system
    TS_ASSERT_THROWS_NOTHING( fts.insert_jump( 5, 15 ) );
    TS_ASSERT_THROWS_NOTHING( fts.insert_jump( 7, 12 ) );

    utility::vector1< core::Size > cycle;
    std::set< core::Size > cycle_set;

    // Cycle is found.
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 1 ) );
    cycle_set.clear();
    std::copy( cycle.begin(), cycle.end(), std::inserter( cycle_set, cycle_set.end() ) );
    TS_ASSERT( cycle_set.find( 5 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 6 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 10 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 14 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 15 ) != cycle_set.end() );

    // A different start point also yields a cycle. (DFS gives same result.)
    TS_ASSERT_THROWS_NOTHING( cycle = fts.cycle( 17 ) );
    cycle_set.clear();
    std::copy( cycle.begin(), cycle.end(), std::inserter( cycle_set, cycle_set.end() ) );
    TS_ASSERT( cycle_set.find( 5 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 6 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 10 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 14 ) != cycle_set.end() );
    TS_ASSERT( cycle_set.find( 15 ) != cycle_set.end() );

    //Unresolvable cycles throw instead of infinite-looping
    utility::vector1< core::Size > bias( 20, 0 );
    bias[ 13 ] = 1;
    TS_ASSERT_THROWS( fts.remove_cycles( bias ), core::environment::EXCN_FTSketchGraph );


    //Resolvable cycle is resolvable
    bias = utility::vector1< core::Size >( 20, 0 );
    bias[ 10 ] = 1;
    bias[ 13 ] = 1;

    TS_ASSERT_THROWS_NOTHING( fts.remove_cycles( bias ) );

    TS_ASSERT( fts.has_cut( 10 ) );
    TS_ASSERT( fts.has_cut( 13 ) );

    // Rendered FT is correct
    core::kinematics::FoldTreeOP ft_out;
    TS_ASSERT_THROWS_NOTHING( ft_out = fts.render() );
    TS_ASSERT_EQUALS( ft_out->num_cutpoint(), 2 );
    TS_ASSERT( ft_out->is_cutpoint( 10 ) );
    TS_ASSERT( ft_out->is_cutpoint( 13 ) );
  }

};
