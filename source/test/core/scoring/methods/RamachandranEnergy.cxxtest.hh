// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/RamachandranEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::RamachandranEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

// Unit headers
#include <core/scoring/methods/RamachandranEnergy.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// Package Headers


#include <numeric/conversions.hh>


// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class RamachandranEnergyTests : public CxxTest::TestSuite {

	public:

	PoseOP the_pose;
	RamachandranEnergyOP rama_energy;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		using namespace std;
		using namespace core;
		core_init();


		the_pose = create_test_in_pdb_poseop();
		//core::import_pose::pose_from_pdb( *the_pose, "core/scoring/methods/test_in.pdb" );

		rama_energy = RamachandranEnergyOP( new RamachandranEnergy );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
		rama_energy.reset();
	}


	// --------------- Test Cases --------------- //
	void dont_test_eval_energy()
	{
		// correct answers taken from rosetta++ v17084
		Real correct_answers[] = { 0,
			1.8437,  -0.246552, 0.0468803,  -0.306571,
			-0.454755,  3.13757, 3.71225,  0.576329,
			-0.435698, -0.942344, 0.255225, -1.41592,
			0.30198, 10.4935, 20, 20, 20, -1.26224, -0.808979,
			20, -0.673049, -0.339168, 0.976442, 0.663541,
			-0.964809, 0.0895034, 13.3146, 2.04172,
			4.49344, 6.30348, -1.28329, -0.513646,
			-0.847239, 0.683858, 2.22122, -0.109381,
			0.830199, 1.12878, -0.0534803, -0.255095,
			-1.19867, -0.76122, -0.356984, -1.29865,
			-0.158319, -1.04536, 0.679841, 3.11838,
			-1.18327, 0.0437514, 0.0696214, -0.985281,
			-0.973587, -0.695756, -1.09403, -0.414759,
			-0.40325, -1.07938, -0.378227, -0.738479,
			-1.0468, -1.05351, -0.714108, -0.0145066,
			-0.286735, -0.0568444, -0.874918, -0.287769,
			-1.1209, -0.296159, -0.564195, -0.654011,
			-0.43274, 0.379993, -1.01306, -1.18156,
			-1.01892, -0.681848, -0.858865, -0.859512,
			1.00974, -0.549346, -0.269065, -0.158142,
			-0.370144, 0.486582, -0.87512, 0.0996326,
			0.160583, -0.828564, -0.310471, -0.0358784,
			0.226202, 0.477853, 0.810734, -0.146505,
			0.876076, 0.0573208, -0.530587, -0.704854,
			-0.974999, 0.905538, 0.247923, 3.98163,
			1.6907, 0.932702, 0.961533, -0.326284,
			0.917562, 0.952441, -1.39607, -0.412053,
			-0.410412, 20, 0 };

	float const TOLERATED_ERROR = 0.001;

		EnergyMap emap;
		for ( int ii = 1; ii <= 116; ++ii ) {
			emap.zero();
			rama_energy->residue_energy( the_pose->residue(ii), *the_pose, emap );
			//std::cout << "residue " << ii << " emap[ rama ] = " << emap[ rama ] << std::endl;
			TS_ASSERT_DELTA( emap[ rama ], correct_answers[ ii - 1 ], TOLERATED_ERROR );
		}


	}


	// --------------- Test Cases --------------- //
	void dont_test_eval_deriv() {

		// correct answers taken from rosetta++ v17084
		Real correct_dE_dphi[] = {
			0, 0.415667, -0.0171281,
			-0.0138161, 0.0410879, 0.0839377, -0.242255,
			0.354805, 0.0703403, 0.0339204, 0.0124663,
			-0.0103839, 0.0503213, 0.11896, -0.663449,
			0, 0, 0, -0.0246639,
			0.0425426, 0, 0.00530617, 0.116217,
			-0.0266832, 0.125175, -0.042308, 0.0161499,
			1.83523, 0.105136, 0.141216, 0.333824,
			-0.0372186, 0.0752653, 0.00821022, -0.0248512,
			0.419967, 0.0235531, 0.0741497, 0.126215,
			0.0443564, -0.032175, 0.0527203, -0.0416654,
			-0.033939, -0.0944868, 0.00102054, 0.0266487,
			0.087497, 0.338022, 0.0168999, -0.0431646,
			0.0496313, -0.0304149, 0.0466007, -0.0417639,
			-0.0179814, -0.0338635, -0.00944268, -0.0322939,
			-0.0187734, -0.0428593, 0.0043185, -0.0400648,
			0.0615923, -0.0353174, -0.0443526, -0.0283583,
			-0.0546665, 0.0822077, 0.0818296, 0.059162,
			-0.0138576, -0.00977368, -0.034723, -0.0092108,
			0.0356659, 0.00155308, 0.0104925, 0.0685211,
			-0.0225038, 0.0411625, -0.172706, -0.00522198,
			0.0511977, 0.0402928, -0.0200223, 0.0588123,
			-0.0352177, 0.0240481, 0.000531677, 0.0471866,
			0.00101169, -0.00992427, -0.00405201, 0.0525659,
			-0.11423, 0.073373, -0.0226876, -0.0245909,
			0.0800682, 0.0420498, -0.0315993, 0.0311541,
			-0.0458183, 0.100521, 0.312332, 0.168934,
			-0.0735353, 0.0835779, 0.111229, 0.0600239,
			-0.0116661, -0.0187792, 0.00339337, 0, 0
		};

		// correct answers taken from rosetta++ v17084
		Real correct_dE_dpsi[] = { 0,
			-0.0460572,  -0.0226938,  -0.0548044, -0.0242637,
			0.0360179, -0.26566, -0.0747737, -0.0669048,
			-0.0424023, 0.0139555, -0.0710771, -0.0142297,
			0.0682576, 0.389362, 0, 0,
			0, 0.0554168, 0.0457162, 0,
			0.0672304, -0.0118136, -0.179038, 0.020087,
			-0.0638011, -0.052335, -0.484653, -0.19253,
			0.269042, 0.247542, -0.0141333, -0.0183849,
			-0.0476692, -0.0467547, 0.0556965, -0.0895534,
			-0.0439458, -0.183755, 0.0022596, -0.0920964,
			-0.000422955, -0.0375434, -0.020991, -0.056979,
			0.00402893, -0.0544251, -0.0566326, -0.174719,
			0.00212482, -0.0641182, -0.0242846, -0.041884,
			-0.0411343, -0.0610899, -0.0458404, -0.0406271,
			-0.0388655, -0.0359497, -0.00117237, -0.0637673,
			-0.0241292, -0.0191541, 0.0125353, -0.0170338,
			0.0158002, 0.0471322, 0.00903028, -0.0789372,
			-0.00208552, 0.003578, -0.027118, -0.0113176,
			-0.0379545, -0.128702, -0.0358132, -0.0198114,
			0.00617447, 0.0387012, 0.00471103, 0.0635429,
			0.116499, -0.0543219, -0.0319358, -0.0389804,
			-0.0406049, -0.0526303, -0.0473781, -0.067742,
			-0.0514289, -0.0265138, -0.0444183, -0.0258162,
			-0.0498829, -0.0379622, -0.0748167, 0.0883636,
			-0.0884887, -0.0411505, -0.0171646, 0.00124844,
			-0.0413568, 0.0568278, -0.0498843, -0.429412,
			-0.148481, -0.0702754, 0.0344602, 0.0743334,
			-0.055903, -0.0526745, 0.0203157, -0.0291913,
			0.0194366, 0, 0
		};

		float const TOLERATED_ERROR = 0.0001;

		ScoreFunction sfxn; // unused in eval deriv
		id::DOF_ID dummy_dofid;
		EnergyMap weights;
		weights[ rama ] = 1.0;
		for ( Size ii = 1; ii <= 116; ++ii ) {
			id::TorsionID ii_phi( ii, id::BB, 1 );
			id::TorsionID ii_psi( ii, id::BB, 2 );
			id::TorsionID ii_omega( ii, id::BB, 3 );

			EnergyDerivative de_dphi = rama_energy->eval_dof_derivative(
				dummy_dofid, ii_phi, *the_pose, sfxn, weights);
			EnergyDerivative correct_de_dphi_in_rad = numeric::conversions::degrees(correct_dE_dphi[ ii - 1 ]);

			Real diff = std::abs( de_dphi - correct_de_dphi_in_rad);
			if ( std::abs( correct_de_dphi_in_rad ) > 1.0 ) {
				diff /= numeric::conversions::degrees( correct_de_dphi_in_rad );
			}
			TS_ASSERT_DELTA( diff, 0, TOLERATED_ERROR );


			EnergyDerivative de_dpsi = rama_energy->eval_dof_derivative(
				dummy_dofid, ii_psi, *the_pose, sfxn, weights);
			EnergyDerivative correct_de_dpsi_in_rad = numeric::conversions::degrees(correct_dE_dpsi[ ii - 1 ]);

			diff = std::abs( de_dpsi - correct_de_dpsi_in_rad);
			if ( std::abs( correct_de_dpsi_in_rad ) > 1.0 ) {
				diff /= numeric::conversions::degrees( correct_de_dpsi_in_rad );
			}
			TS_ASSERT_DELTA( diff, 0, TOLERATED_ERROR );

			EnergyDerivative de_domega = rama_energy->eval_dof_derivative(
				dummy_dofid, ii_omega, *the_pose, sfxn, weights);
			TS_ASSERT_DELTA( de_domega, 0, TOLERATED_ERROR );

		}

	}


	void test_setup_for_minimizing_with_rama_and_full_bb_flex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;


		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama,   0.4 );

		kinematics::MoveMap movemap; movemap.set_bb( 10, true );movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}

	void dont_test_atom_tree_minimize_with_rama_and_etable()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( rama,   0.4 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -27.89673840790607, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -28.10526841844615, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_setup_for_minimizing_with_rama_and_partial_bbflex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama,   0.4 );

		kinematics::MoveMap movemap; movemap.set_bb( 10, true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}

	void dont_test_atom_tree_minimize_with_rama_and_etable2()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( rama,   0.4 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -27.89673840790607, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -27.9539686076814, end_score, 1e-12 );
	}

	void test_rama_deriv_check()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama,   0.4 );
		kinematics::MoveMap movemap; movemap.set_bb( true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void dont_test_write_rama_energies() {
		using namespace core::scoring;
		Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
		// lets output a table of values from the leucine phi/psi distribution
		for ( Size ii = 0; ii < 360; ++ii ) {
			Real ii_phi = -180.0 + ii;
			for ( Size jj = 0; jj < 360; ++jj ) {
				Real jj_psi = -180.0 + jj;
				std::cout << rama.eval_rama_score_residue( core::chemical::aa_leu, ii_phi, jj_psi );
				if ( jj != 359 ) std::cout << " ";
			}
			std::cout << "\n";
		}
	}

};


