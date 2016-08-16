// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/GrowRight.cxxtest.hh
/// @brief  unit tests for GrowRight BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/GrowRight.hh>

#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class GrowRightTests : public CxxTest::TestSuite
{


private: // data


	core::pose::Pose pose_;


public: // setup


	typedef std::string String;
	typedef protocols::forge::build::GrowRight GrowRight;
	typedef protocols::forge::build::Interval Interval;


	GrowRightTests() {};


	// Shared initialization.
	void setUp() {
		core_init();

		// create dummy pose
		core::pose::make_pose_from_sequence(
			pose_,
			"A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose_.n_residue(); i <= ie; ++i ) {
			pose_.set_secstruct( i, 'L' );
		}
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests


	/// @brief test C-terminal extension
	void test_c_term_extension() {
		core::Size ori_n_residue = pose_.n_residue();

		GrowRight grow( pose_.n_residue(), String( 7, 'H' ), "AAAAAAG" ); // 7-mer, full-atom
		grow.modify( pose_ );

		TS_ASSERT_EQUALS( grow.pos(), ori_n_residue );
		TS_ASSERT_EQUALS( grow.original_interval().left, ori_n_residue );
		TS_ASSERT_EQUALS( grow.original_interval().right, ori_n_residue );
		TS_ASSERT( !grow.original_interval_valid() );
		TS_ASSERT_EQUALS( pose_.n_residue(), 27 );
		TS_ASSERT_EQUALS( pose_.fold_tree().num_cutpoint(), 0 );
		TS_ASSERT( pose_.residue( pose_.n_residue() ).is_upper_terminus() );
		TS_ASSERT_EQUALS( pose_.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWYAAAAAAG[GLY:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose_.secstruct(), "LLLLLLLLLLLLLLLLLLLLHHHHHHH" );
	}


	/// @brief test internal C-side extension
	void test_c_side_internal() {
		GrowRight grow( 9, String( 7, 'H' ), "AAAAAAG" ); // 7-mer, full-atom
		grow.modify( pose_ );

		TS_ASSERT_EQUALS( grow.pos(), 9 );
		TS_ASSERT_EQUALS( grow.original_interval().left, 9 );
		TS_ASSERT_EQUALS( grow.original_interval().right, 9 );
		TS_ASSERT( !grow.original_interval_valid() );
		TS_ASSERT_EQUALS( pose_.n_residue(), 27 );
		TS_ASSERT_EQUALS( pose_.fold_tree().num_cutpoint(), 0 );
		TS_ASSERT_EQUALS( pose_.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKAAAAAAGLMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose_.secstruct(), "LLLLLLLLLHHHHHHHLLLLLLLLLLL" );
	}


};
