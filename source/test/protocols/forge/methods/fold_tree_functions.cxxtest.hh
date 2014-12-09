// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/forge/build/GrowLeft.cxxtest.hh
/// @brief  unit tests for GrowLeft BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>

//Auto Headers
#include <utility/vector1.hh>


// non-camel case to be consistent with library function file
class fold_tree_functions_Tests : public CxxTest::TestSuite
{


public: // setup


	typedef core::Size Size;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;


	fold_tree_functions_Tests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests


	/// @brief test N-terminal extension
	/// @remarks We only need to check GrowLeft here, as NtermExt changes
	///  none of the engine.
	void test_remove_cutpoints() {
		using protocols::forge::methods::remove_cutpoints;

		// fake a fold tree with 3 cutpoints
		FoldTree ft;
		ft.simple_tree( 50 );
		ft.new_jump( 3, 10, 7 );
		ft.new_jump( 14, 26, 20 );
		ft.new_jump( 35, 47, 41 );
		ft.reorder( 1 );

		// remove cutpoints at 7 and 41
		utility::vector1< Size > cutpoints;
		cutpoints.push_back( 7 );
		cutpoints.push_back( 41 );

		remove_cutpoints( cutpoints, ft );
		TS_ASSERT( !ft.is_cutpoint( 7 ) );
		TS_ASSERT( ft.is_cutpoint( 20 ) );
		TS_ASSERT( !ft.is_cutpoint( 41 ) );
	}



};
