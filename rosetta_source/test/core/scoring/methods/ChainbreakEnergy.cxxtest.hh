// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/ChainbreakEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::ChainbreakEnergy
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/ChainbreakEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class ChainbreakEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	/// This can be used to ensure norm matches norm-numeric
	void test_chainbreak_deriv_check()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose = pdb1ubq5to13_pose();
		FoldTree ft;
		ft.add_edge( 1, 3, -1 );
		ft.add_edge( 3, 7,  1 );
		ft.add_edge( 7, 9, -1 );
		ft.add_edge( 7, 6, -1 );
		ft.add_edge( 3, 5, -1 );

		pose.fold_tree( ft );
		pose.set_phi( 4, 180 );
		pose.set_psi( 4, 180 );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, 5 );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, 6 );


		//pose.dump_pdb( "test_ft1.pdb" );

		ScoreFunction sfxn;
		sfxn.set_weight( chainbreak, 0.5 );
		TS_ASSERT_DELTA( sfxn( pose ), 126.855, 1e-3 );

		//std::cout << "chainbreak score: " << sfxn( pose ) << std::endl;
		MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

};


