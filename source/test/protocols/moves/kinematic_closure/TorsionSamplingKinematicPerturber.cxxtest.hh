// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/kinematic_closure/TorsionSamplingKinematicPerturber.cxxtest.hh
/// @brief  test for TorsionSamplingKinematicPerturber
/// @author Steven Lewis

//This is really a test for the TorsionSamplingKinematicPerturber, but for a variety of reasons it is vastly easier to pass it through the KinematicMover than test that class directly.  For the rama-sampling and local-sampling modes, it tests that the output poses' torsions match what they did when the test was created, and for the local sampling mode that the results are indeed local to the starting point.

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

class TorsionSamplingKinematicPerturberTests : public CxxTest::TestSuite {

	core::pose::Pose pose;
	core::Size start, end, middle, nres;
	utility::vector1<core::Real> torsions;

	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kin_mover;
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP kin_mover_cap;
	protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP TsamplingKP;
	/* AS Oct 04, 2012 -- note that there currently is no test for the actual TorsionRestrictedKinematicPerturber, due to OS-dependent instability issues */

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	TorsionSamplingKinematicPerturberTests(){
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

		TsamplingKP = protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP( new protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber(kin_mover_cap) );
		kin_mover->set_perturber(TsamplingKP);
	}

	virtual ~TorsionSamplingKinematicPerturberTests() {}

	static TorsionSamplingKinematicPerturberTests* createSuite() {
		return new TorsionSamplingKinematicPerturberTests();
	}

	static void destroySuite( TorsionSamplingKinematicPerturberTests *suite ) {
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
	//This test deactivated; unreliable on mac versus linux.  You can use the extra-commented-out Ramachandran code to notice that the Ramanchandran object random_phipsi_from_rama is getting called different numbers of times on the two architectures; this is symptomatic of KIC taking different trajectories.  This test is disabled until someone with access to both architectures chooses to debug it further.
	/*	void te //line broken to fix cxxtest detection
			st_TorsionSamplingKinematicPerturber_rama(){

		//set up kinematic mover
		kin_mover->set_vary_bondangles( false );
		kin_mover->set_sample_nonpivot_torsions( true );
		kin_mover->set_rama_check( false );
		kin_mover->set_idealize_loop_first( false );
		kin_mover->set_pivots(start, middle, end);

		//set up kinematic perturber
		TsamplingKP->set_vary_ca_bond_angles(false);
		TsamplingKP->set_sample_vicinity(false);

		core::pose::Pose copy(pose);
		//copy.dump_pdb("1.pdb");

		kin_mover->apply(copy);

		//		std::cout << std::boolalpha << kin_mover->last_move_succeeded() << std::endl;

		//copy.dump_pdb("2.pdb");

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

		//This part of the protein is different; here's raw torsions
		//7
		TS_ASSERT_DELTA( 287.38054939225827411065, copy.phi(7), 1e-9 );
		TS_ASSERT_DELTA( 342.29681462866511765242, copy.psi(7), 1e-9 );
		TS_ASSERT_DELTA( -176.72579983092379052323, copy.omega(7), 1e-9 );
		//8
		TS_ASSERT_DELTA( -59.21399219966150440086, copy.phi(8), 1e-9 );
		TS_ASSERT_DELTA( -29.34088040118307105786, copy.psi(8), 1e-9 );
		TS_ASSERT_DELTA( 178.13567298531563665165, copy.omega(8), 1e-9 );
		//9
		TS_ASSERT_DELTA( 256.48375207684176757539, copy.phi(9), 1e-9 );
		TS_ASSERT_DELTA( 36.93239779579432280343, copy.psi(9), 1e-9 );
		TS_ASSERT_DELTA( 173.15615572009770062323, copy.omega(9), 1e-9 );
		//10
		TS_ASSERT_DELTA( 77.03681325271357138718, copy.phi(10), 1e-9 );
		TS_ASSERT_DELTA( -10.85342736709987221388, copy.psi(10), 1e-9 );
		TS_ASSERT_DELTA( 175.79932009520385349788, copy.omega(10), 1e-9 );
		//11
		TS_ASSERT_DELTA( 89.01820885863949683880, copy.phi(11), 1e-9 );
		TS_ASSERT_DELTA( 252.76124586141020245122, copy.psi(11), 1e-9 );
		TS_ASSERT_DELTA( -178.41640975783539602162, copy.omega(11), 1e-9 );

		//This part of the protein is post-loop and should not change
		check_similar_parts(end+1, nres, copy);

// 		//hack in random_phipsi_from_rama test

// 		core::scoring::Ramachandran rama_obj;
// 		core::Real phi, psi;

// 		//print code comparison lines
// 		// 		std::cout << std::setprecision(20) << std::fixed;
// 		// 		for(int i(1); i<15; ++i){
// 		// 			rama_obj.random_phipsi_from_rama(core::chemical::aa_ala, phi, psi);
// 		// 			std::cout << "rama_obj.random_phipsi_from_rama(core::chemical::aa_ala, phi, psi);" << std::endl;
// 		// 			std::cout << "TS_ASSERT_DELTA( " << phi << ", phi, 1e-9 );" << std::endl;
// 		// 			std::cout << "TS_ASSERT_DELTA( " << psi << ", psi, 1e-9 );" << std::endl << std::endl;
// 		// 		}

// 		std::cout << "Running the new extra tests" << std::endl;

// 		//comparisons go here

		return;
	}*/

	void test_TorsionSamplingKinematicPerturber() {
		
		// this test is empty... the compiler complains if there's no test defined
		return;
	}

};//end class
