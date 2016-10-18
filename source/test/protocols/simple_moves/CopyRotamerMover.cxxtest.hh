// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/membrane/CopyRotamerMover.cxxtest.hh
/// @brief   Unit test for the CopyRotamerMover.
/// @author  Vikram K. Mulligan (vmullig@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// Package Headers
#include <protocols/simple_moves/CopyRotamerMover.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer TR("protocols.simple_moves.CopyRotamerMover.cxxtest");

class CopyRotamerMoverTests : public CxxTest::TestSuite {

private: // variables
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {

		using namespace core::pose;
		using namespace core::import_pose;

		// Initialize core & options system
		core_init();

		// Load in pose from pdb
		testpose_ = pdb1ubq5to13_poseop();
		core::scoring::ScoreFunctionOP scorefxn_;

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////

	/// @brief Test the CopyRotamerMover
	///
	void test_copyrotamermover() {
		using namespace protocols::simple_moves;
		core::pose::PoseOP testpose_copy_(testpose_->clone());
		//testpose_copy_->dump_pdb("before.pdb");//DELETE ME

		utility::vector1 <core::Real> res4_chivals;
		res4_chivals = testpose_->residue(4).chi();

		//Directly mutate L8F
		CopyRotamerMoverOP copyrot( new CopyRotamerMover );
		copyrot->set_template_res_index(4);
		copyrot->set_target_res_index(8);
		copyrot->apply(*testpose_copy_);


		TS_ASSERT_EQUALS( res4_chivals.size() , 2  );
		TS_ASSERT_EQUALS( testpose_copy_->residue(8).nchi() , 2  );

		TS_ASSERT_EQUALS( testpose_copy_->residue_type(4).name(), testpose_copy_->residue_type(8).name() );

		TS_ASSERT_DELTA( testpose_copy_->chi(1,4), testpose_copy_->chi(1,8), 0.0001  );
		TS_ASSERT_DELTA( testpose_copy_->chi(2,4), testpose_copy_->chi(2,8), 0.0001  );

		return;
	}


}; // MutateResidueTests
