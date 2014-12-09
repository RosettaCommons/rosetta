// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/kinematic_closure/VicinitySamplingKinematicPerturber.cxxtest.hh
/// @brief  test for VicinitySamplingKinematicPerturber
/// @author Steven Lewis
/// @author Amelie Stein amelie.stein@ucsf.edu (implementation update)

//This is really a test for the VicinitySamplingKinematicPerturber, but for a variety of reasons it is vastly easier to pass it through the KinematicMover than test that class directly.  It tests that the output poses' torsions match what they did when the test was created, and that the results are indeed local to the starting point.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/pose_funcs.hh>

// Unit header
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
// AUTO-REMOVED #include <protocols/loops/KinematicWrapper.hh>

// project headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>

// #include <core/scoring/Ramachandran.hh>
// #include <core/chemical/AA.hh>

// AUTO-REMOVED #include <numeric/xyz.io.hh>

//Auto Headers
#include <utility/vector1.hh>

// --------------- Test Class --------------- //

class VicinitySamplingKinematicPerturberTests : public CxxTest::TestSuite {

	core::pose::Pose pose;
	core::Size start, end, middle, nres;
	utility::vector1<core::Real> torsions;

	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kin_mover;
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP kin_mover_cap;
	protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP VsamplingKP;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	VicinitySamplingKinematicPerturberTests(){
		core_init();
		pose = create_trpcage_ideal_pose();
		//pose.dump_pdb("initialpdb.pdb");

		//I checked; these work on the first try
		start=7; end=11; middle=9; nres=pose.total_residue();
		torsions.reserve(nres*3);

		//extract original torsions
		for(core::Size i(1); i<=nres; ++i){
			torsions.push_back(pose.phi(i));
			torsions.push_back(pose.psi(i));
 			torsions.push_back(pose.omega(i));
		}

		kin_mover = protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP( new protocols::loops::loop_closure::kinematic_closure::KinematicMover() );
		kin_mover_cap = protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP( kin_mover );

		VsamplingKP = protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP( new protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber(kin_mover_cap) );
		kin_mover->set_perturber(VsamplingKP);
	}

	virtual ~VicinitySamplingKinematicPerturberTests() {}

	static VicinitySamplingKinematicPerturberTests* createSuite() {
		return new VicinitySamplingKinematicPerturberTests();
	}

	static void destroySuite( VicinitySamplingKinematicPerturberTests *suite ) {
		delete suite;
	}

	void setUp() {}

	void tearDown() {}



	// ------------- Helper Functions ------------- //

	void check_similar_parts (core::Size s, core::Size e, core::pose::Pose copy){
		for(core::Size i(s); i<=e; ++i){
			//std::cout << i << std::endl;
			TS_ASSERT_DELTA( torsions[(i*3)-2], copy.phi(i), 1e-9 );
			TS_ASSERT_DELTA( torsions[(i*3)-1], copy.psi(i), 1e-9 );
			TS_ASSERT_DELTA( torsions[(i*3)-0], copy.omega(i), 1e-9 );
		}
	}

	// --------------- Test Cases --------------- //

	void test_VicinitySamplingKinematicPerturber(){

		//set up kinematic mover
		kin_mover->set_vary_bondangles( false );
		kin_mover->set_sample_nonpivot_torsions( true );
		kin_mover->set_rama_check( false );
		kin_mover->set_idealize_loop_first( false );
		kin_mover->set_pivots(start, middle, end);

		//set up kinematic perturber
		VsamplingKP->set_vary_ca_bond_angles(false);

		core::pose::Pose copy(pose);
		//copy.dump_pdb("1a.pdb");

		kin_mover->apply(copy);

		//copy.dump_pdb("2a.pdb");

		//This part of the protein is pre-loop and should not change
		check_similar_parts(1, start-1, copy);

		//This code will print the torsions for the loop
// 		std::cout << std::setprecision(20) << std::fixed;
// 		for(core::Size i(start); i<=end; ++i){
// 			std::cout << "//" << i << std::endl;
// 			std::cout << "TS_ASSERT_DELTA( " << copy.phi(i) << ", copy.phi(" << i << "), 1e-9 );\n";
// 			std::cout << "TS_ASSERT_DELTA( " << copy.psi(i) << ", copy.psi(" << i << "), 1e-9 );\n";
// 			std::cout << "TS_ASSERT_DELTA( " << copy.omega(i) << ", copy.omega(" << i << "), 1e-9 );\n";
// 		}

		test::UTracer UT("protocols/moves/kinematic_closure/VicinitySamplingKinematicPerturber.u");
		UT.abs_tolerance(1e-9);

		UT << "\nRaw torsions:\n";
		for(int i=7; i<12; i++) {
			UT << "  copy.phi(" << i << "):   " << copy.phi(i) << "\n";
			UT << "  copy.psi(" << i << "):   " << copy.psi(i) << "\n";
			UT << "  copy.omega(" << i << "): " << copy.omega(i)<< "\n";
		}

		/*
		//This part of the protein is different; here's raw torsions
		//7
		TS_ASSERT_DELTA( 289.40271436608145450009, copy.phi(7), 1e-9 );
		TS_ASSERT_DELTA( 318.44772587254220752584, copy.psi(7), 1e-9 );
		TS_ASSERT_DELTA( -176.72579983092379052323, copy.omega(7), 1e-9 );
		//8
		TS_ASSERT_DELTA( -62.57695096442053284136, copy.phi(8), 1e-9 );
		TS_ASSERT_DELTA( -18.86871213525618173890, copy.psi(8), 1e-9 );
		TS_ASSERT_DELTA( 178.13567298531563665165, copy.omega(8), 1e-9 );
		//9
		TS_ASSERT_DELTA( 272.60910140670944201702, copy.phi(9), 1e-9 );
		TS_ASSERT_DELTA( 352.05672034257025870829, copy.psi(9), 1e-9 );
		TS_ASSERT_DELTA( 173.15615572009770062323, copy.omega(9), 1e-9 );
		//10
		TS_ASSERT_DELTA( 113.66515028602398729163, copy.phi(10), 1e-9 );
		TS_ASSERT_DELTA( 11.94382251408127437742, copy.psi(10), 1e-9 );
		TS_ASSERT_DELTA( 175.79932009520385349788, copy.omega(10), 1e-9 );
		//11
		TS_ASSERT_DELTA( 58.84742871607191005978, copy.phi(11), 1e-9 );
		TS_ASSERT_DELTA( 234.70105204067374415899, copy.psi(11), 1e-9 );
		TS_ASSERT_DELTA( -178.41640975783539602162, copy.omega(11), 1e-9 );
		*/

		// With "local sampling", it should be within one degree of where it started, except
		// that the perturbation is scaled by a Gaussian random number, so it could be +/-
		// 180 degrees.  *sigh*
		for(core::Size i(start+1); i<=end-1; ++i) {
			if( i==middle) continue; //this is a pivot; not handled by TorsionSamplingKinematicPerturber
			TS_ASSERT_DELTA( torsions[(i*3)-2], copy.phi(i), 2 ); //delta to two degrees!
			TS_ASSERT_DELTA( torsions[(i*3)-1], copy.psi(i), 2 );
			TS_ASSERT_DELTA( torsions[(i*3)-0], copy.omega(i), 2 );
		}

		//This part of the protein is post-loop and should not change
		check_similar_parts(end+1, nres, copy);

		return;
	}

};//end class
