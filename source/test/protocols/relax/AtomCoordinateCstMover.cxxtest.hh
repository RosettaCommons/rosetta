// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/relax/AtomCoordinateCstMover.cxxtest.hh
/// @brief  test for AtomCoordinateCstMover
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pdb1rpb.hh>

// Unit header
#include <protocols/relax/AtomCoordinateCstMover.hh>

//#include <protocols/relax/FastRelax.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

// project headers
#include <core/types.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//class OldAddCst : public protocols::relax::FastRelax {
//public:
//	OldAddCst() {
//		set_scorefxn( core::scoring::getScoreFunction() );
//	}
//	void apply( core::pose::Pose & pose ) {
//		core::kinematics::MoveMapOP local_movemap = get_movemap()->clone();
//		initialize_movemap( pose, *local_movemap );
//		set_movemap(local_movemap);
//		set_up_constraints( pose, *local_movemap );
//	}
//};

 // --------------- Test Class --------------- //

class AtomCoordinateCstMoverTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	void test_input_harmonic_default() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose(pdb1rpb_pose());
		core::pose::addVirtualResAsRoot(pose);

		//option[ OptionKeys::relax::constrain_relax_to_start_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( false );
		//option[ OptionKeys::relax::coord_cst_stdev ].value( 0.5 );
		//option[ OptionKeys::relax::coord_cst_width ].to_default();
		//OldAddCst Cst_mover;

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.5);

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_input_harmonic_default.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

	void test_input_bounded_default() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose(pdb1rpb_pose());
		core::pose::addVirtualResAsRoot(pose);

		//option[ OptionKeys::relax::constrain_relax_to_start_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( false );
		//option[ OptionKeys::relax::coord_cst_stdev ].value( 0.5 );
		//option[ OptionKeys::relax::coord_cst_width ].value( 0.25 );
		//OldAddCst Cst_mover;

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.5);
		Cst_mover.bounded( true );
		Cst_mover.cst_width( 0.25 );

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_input_bounded_default.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

	void test_native_harmonic_default() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose(pdb1rpb_pose());
		core::pose::addVirtualResAsRoot(pose);

		core::pose::PoseOP native = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native, "protocols/relax/AtomCoordinateCstMover_native.pdb" );

		//option[ OptionKeys::relax::constrain_relax_to_native_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( false );
		//option[ OptionKeys::relax::coord_cst_stdev ].value( 0.5 );
		//option[ OptionKeys::relax::coord_cst_width ].to_default();
		//OldAddCst Cst_mover;
		//Cst_mover.set_native_pose(native);

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.5);
		Cst_mover.set_refstruct( native );

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_native_harmonic_default.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

	void test_input_harmonic_sidechain() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose(pdb1rpb_pose());
		core::pose::addVirtualResAsRoot(pose);

		//option[ OptionKeys::relax::constrain_relax_to_start_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( true );
		//option[ OptionKeys::relax::coord_cst_stdev ].value( 0.75 );
		//option[ OptionKeys::relax::coord_cst_width ].to_default();
		//OldAddCst Cst_mover;

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.75);
		Cst_mover.cst_sidechain( true );

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_input_harmonic_sidechain.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

	void test_native_harmonic_sidechain() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose(pdb1rpb_pose());
		core::pose::addVirtualResAsRoot(pose);

		core::pose::PoseOP native = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native, "protocols/relax/AtomCoordinateCstMover_native.pdb" );

		//option[ OptionKeys::relax::constrain_relax_to_native_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( true );
		//option[ OptionKeys::relax::coord_cst_stdev ].value( 0.75 );
		//option[ OptionKeys::relax::coord_cst_width ].to_default();
		//OldAddCst Cst_mover;
		//Cst_mover.set_native_pose(native);

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.75);
		Cst_mover.cst_sidechain( true );
		Cst_mover.set_refstruct( native );

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_native_harmonic_sidechain.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

	void test_integration_native() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::relax;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/relax/1a19.pdb" );

		core::pose::addVirtualResAsRoot(pose);

		core::pose::PoseOP native = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native, "protocols/relax/1a19_trunc.pdb" );

		//option[ OptionKeys::relax::constrain_relax_to_native_coords ].value( true );
		//option[ OptionKeys::relax::coord_constrain_sidechains ].value( true );
		//option[ OptionKeys::relax::coord_cst_stdev ].to_default();
		//option[ OptionKeys::relax::coord_cst_width ].to_default();
		//OldAddCst Cst_mover;
		//Cst_mover.set_native_pose(native);

		AtomCoordinateCstMover Cst_mover;
		Cst_mover.cst_sd(0.5);
		Cst_mover.cst_sidechain( true );
		Cst_mover.set_refstruct( native );

		Cst_mover.apply(pose);

		test::UTracer UT("protocols/relax/AtomCoordinateCstMover_integration_native.cst");
		pose.constraint_set()->show_definition(UT,pose);
	}

};//end class
