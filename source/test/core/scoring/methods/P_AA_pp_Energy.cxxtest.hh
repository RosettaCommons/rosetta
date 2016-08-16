// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/RamachandranEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::RamachandranEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/P_AA_pp_Energy.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/ScoringManager.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>


#include <numeric/conversions.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class P_AA_pp_EnergyTests : public CxxTest::TestSuite {

public:

	PoseOP the_pose;
	P_AA_pp_EnergyOP paapp_energy;


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
		//core::import_pose::pose_from_file( *the_pose, "core/scoring/methods/test_in.pdb" , core::import_pose::PDB_file);

		paapp_energy = P_AA_pp_EnergyOP( new P_AA_pp_Energy );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
		paapp_energy.reset();
	}


	// --------------- Test Cases --------------- //
	void dont_test_eval_energy()
	{
		// correct answers taken from rosetta++ v17084
		Real correct_answers[] = { 0,
			-0.267003, -0.172068, -0.0234335, -0.0781696, 0.183428,
			0.393806, 0.132634, -0.272525, 0.160525, 0.0693895,
			-0.470709, -0.409499, -0.990583, -0.778315, -1.15476,
			-0.888028, -1.24405, -2.03641, -0.114451, -1.17934,
			0.193517, -0.34491, -0.472938, 0.445191, -0.0738133,
			-0.325311, -0.313286, 0.704141, -2.4589, 0.529777,
			-0.422739, -0.149188, 0.0406586, -0.0182416, 6.66998e-05,
			-0.152545, -0.272386, 0.0610415, 0.230158, -0.138024,
			-0.123104, -0.499217, -1.14451, -2.04541, 0.150509,
			0.153041, -0.25365, 0.699729, -0.185124, 0.171611,
			-0.179547, -0.26841, 0.0590553, -0.267694, -0.125411,
			-0.282155, 0.365959, -0.46059, 0.34264, -0.299417,
			-0.362812, -0.235615, -0.0747837, -0.231082, 0.0652427,
			-1.03201, -0.827206, -0.103322, -2.0602, 0.4331,
			-0.367151, 0.346117, -0.33369, -0.166156, -0.348641,
			-0.279444, -0.374949, -2.52621, -0.231412, -0.313947,
			-0.563522, -0.283233, -0.175948, 0.332023, -0.343283,
			0.320543, -0.253355, 0.249818, -0.307421, -0.139671,
			0.177677, -0.0319318, -0.250677, 0.609964, -0.0434868,
			-1.91349, 1.40635, -0.00860932, -0.156088, 0.398539,
			-0.252835, -0.133606, -0.273483, 0.123033, 0.353364,
			0.0908511, -1.81416, -0.362164, -0.0563304, 1.02746,
			0.105663, -0.413347, 0.0460336, -2.80756,
			0
			};

		float const TOLERATED_ERROR = 0.0001;

		EnergyMap emap;
		for ( int ii = 1; ii <= 116; ++ii ) {
			emap.zero();
			paapp_energy->residue_energy( the_pose->residue(ii), *the_pose, emap );
			//std::cout << "residue " << ii << " emap[ p_aa_pp ] = " << emap[ p_aa_pp ] << std::endl;
			TS_ASSERT_DELTA( emap[ p_aa_pp ], correct_answers[ ii - 1 ], TOLERATED_ERROR );
		}


	}


	// --------------- Test Cases --------------- //
	void dont_test_eval_deriv() {

		// correct answers taken from rosetta++ v17084
		Real correct_dE_dphi[] = { 0,
			0.00233691, 0.00795944, 0.0392169, 0.0159329, 0.00295341,
			0.00363082, -0.0510909, -0.000245607, -0.00274645, -0.00945298,
			0.0174122, 0.00844788, 0.0137428, -0.0480851, -0.105097,
			0.154839, 0.0568625, -0.0353066, 0.0334361, 0.0727225,
			0.0171216, -0.00894684, -0.000141731, -0.00411311, 0.00903221,
			-0.000925144, -0.00786996, 0.0161266, -0.0023924, 0.107322,
			-0.0172554, -0.0124634, -0.0114224, 0.0118325, -0.00942108,
			0.00423673, -0.000513296, -0.00582425, 0.000590451, 0.013281,
			-0.00205015, 0.00120534, 0.0216519, -0.0630383, 0.000209275,
			-0.00452167, 0.00337774, -0.0434146, 0.00648143, -0.00194848,
			0.0267815, 0.00233581, -0.00147893, 0.00188727, 0.00208954,
			-0.00123905, 0.0153526, 0.0122959, 0.000456536, -0.0386242,
			0.00369934, 0.0254753, 0.0288034, 0.00941142, -0.0329291,
			0.0217585, -0.0109648, 0.00402428, 0.0177926, 0.0130636,
			-0.00228076, 0.0147811, -0.0143062, -0.0121494, 0.0101027,
			0.00762261, 0.0413676, -0.000932433, 0.00247019, -0.0135157,
			-0.0680968, 0.00181827, 0.0137013, 0.0118509, -0.002108,
			-0.00126184, 0.0017647, 0.00855062, -0.00460128, -0.00279352,
			0.00397245, 0.0184209, 0.0194357, 0.179625, -0.0711419,
			-0.0184867, -0.208629, 0.018017, -0.00880259, 0.00441168,
			0.0019097, -0.0110793, -0.0315393, 0.0776582, 0.0687826,
			-0.0278022, -0.0107653, -0.001359, -0.0417697, 0.0982451,
			0.0251439, -0.0346583, 0.00467711, -0.111108, 0
			};

		// correct answers taken from rosetta++ v17084
		Real correct_dE_dpsi[] = { 0,
			0.00279095, 0.049242, 0.00451825, 0.0528715, -0.0136721,
			-0.0683584, -0.0276967, 0.00892394, -0.0161514, -0.00824397,
			0.0158715, 0.00181597, 0.0614149, 0.0483811, 0.0130659,
			0.143131, -0.041473, 0.0127209, 0.0292022, 0.148331,
			0.0643031, -0.00410481, 0.0302769, -0.00753456, -0.0367733,
			-0.0103398, -0.103985, -0.0181595, 0.0577008, 0.0982017,
			-0.0427146, -0.00116093, 0.0187129, -0.0458294, 0.00388936,
			-0.0165709, 0.012839, -0.0365556, 0.0588186, 0.00129469,
			-0.000388775, -0.00777659, 0.0247349, -0.00462179, 0.00860937,
			-0.0239551, 0.0158138, 0.0165468, 0.0144646, 0.0387009,
			0.0343739, -0.00905954, -0.0504798, -0.0111041, 0.000384891,
			-0.0196523, 0.0198873, -0.00864702, 0.0262759, -0.00884102,
			0.0199125, 0.00655608, 0.0241419, 0.0295324, 0.000355545,
			-0.0126152, -0.038484, -0.0148608, -0.0305592, 0.0268758,
			-0.0039788, 0.00188305, -0.00132879, -0.0352135, -0.00631891,
			0.00317642, -0.00980082, 0.00108159, -0.00955604, -0.00025777,
			0.0337909, -0.0186785, -0.0070708, 0.0189615, -0.00352758,
			0.00670957, -0.00949354, 0.000433131, 0.0476124, -0.00179075,
			-0.0213258, 0.0207525, -0.00227411, 0.020654, -0.0295131,
			0.00810109, -0.0424401, 0.0225212, -0.000615853, 0.0042428,
			-0.00869106, -0.0391215, -0.017574, 0.0149176, 0.0559156,
			0.0652406, -0.0272102, 0.0237519, -0.0301465, 0.027841,
			-0.00126092, -0.000286411, 0.0139144, -0.115895, 0
			};

		float const TOLERATED_ERROR = 0.0001;

		ScoreFunction sfxn; // unused in eval deriv
		id::DOF_ID dummy_dofid;
		EnergyMap weights;
		weights[ p_aa_pp ] = 1.0;
		for ( Size ii = 1; ii <= 116; ++ii ) {
			id::TorsionID ii_phi( ii, id::BB, 1 );
			id::TorsionID ii_psi( ii, id::BB, 2 );
			id::TorsionID ii_omega( ii, id::BB, 3 );

			EnergyDerivative de_dphi = paapp_energy->eval_dof_derivative(
				dummy_dofid, ii_phi, *the_pose, sfxn, weights);
			EnergyDerivative correct_de_dphi_in_rad = numeric::conversions::degrees(correct_dE_dphi[ ii - 1 ]);

			Real diff = std::abs( de_dphi - correct_de_dphi_in_rad);
			if ( std::abs( correct_de_dphi_in_rad ) > 1.0 ) {
				diff /= numeric::conversions::degrees( correct_de_dphi_in_rad );
			}
			TS_ASSERT_DELTA( diff, 0, TOLERATED_ERROR );


			EnergyDerivative de_dpsi = paapp_energy->eval_dof_derivative(
				dummy_dofid, ii_psi, *the_pose, sfxn, weights);
			EnergyDerivative correct_de_dpsi_in_rad = numeric::conversions::degrees(correct_dE_dpsi[ ii - 1 ]);

			diff = std::abs( de_dpsi - correct_de_dpsi_in_rad);
			if ( std::abs( correct_de_dpsi_in_rad ) > 1.0 ) {
				diff /= numeric::conversions::degrees( correct_de_dpsi_in_rad );
			}
			TS_ASSERT_DELTA( diff, 0, TOLERATED_ERROR );

			EnergyDerivative de_domega = paapp_energy->eval_dof_derivative(
				dummy_dofid, ii_omega, *the_pose, sfxn, weights);
			TS_ASSERT_DELTA( de_domega, 0, TOLERATED_ERROR );

		}

	}

	void dont_test_paappE_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( p_aa_pp, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -5.845401447356349, false /* print start score?*/ );

	}

	void dont_test_paappE_derivs_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( p_aa_pp, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nonzero_deriv_only( true );

		//adv.compute_pose_atom_derivs();

		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 6, 1), -0.0810506479586785,-0.1729286237967581,0.145681665969538,-0.01926002587211747,0.01822679333460214,0.01092036324788663));
		adl.add( AtomDeriv( AtomID( 7, 5), 0.1451026311926192,0.2751539263066936,-0.2486368093016149,0.04121673430575654,-0.02827994486596643,-0.007242215156046175));
		adl.add( AtomDeriv( AtomID( 9, 8), -0.04355231855740711,-0.2492759577786394,0.2652270496793911,-0.05292169880976719,-0.005537079151122526,-0.01389422156717781));
		adl.add( AtomDeriv( AtomID( 6, 9), -0.06298998034744646,-1.079769123857607,1.164673190929659,-0.2320652225216741,-0.02712997805317142,-0.03770319157399334));
		adl.add( AtomDeriv( AtomID( 6, 13), 0.3775477578615385,0.3165614471106651,-0.04924384360965836,0.01898061306255181,-0.01315886435595045,0.06093144933070874));
		adl.add( AtomDeriv( AtomID( 6, 14), -0.7677617536082132,-0.1564785265613581,-1.087658990121733,0.1192971034081412,0.09145476220290806,-0.09736733729171168));
		adl.add( AtomDeriv( AtomID( 9, 16), 0.4327043114175876,1.066736858577004,-0.1900422635455818,0.1247524964271092,-0.03557568911129938,0.08435515301033364));

		adv.validate_atom_deriv_list( adl );

	}

	void dont_test_paappE_start_score_start_func_match_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( p_aa_pp, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -5.845401447356349,  false /* print start score?*/ );

	}

	void dont_test_paappE_derivs_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( p_aa_pp, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.set_nonzero_deriv_only( true );

		//adv.compute_pose_atom_derivs();

		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 6, 9), -0.04249031567091288,-1.226819779126311,1.326945097276973,-0.2630302128978022,-0.04272020873565824,-0.04791926504933069));
		adl.add( AtomDeriv( AtomID( 6, 14), 0.4590761069158174,0.607139168647146,-0.601577773978451,0.1242259490666155,-0.0236160089768855,0.07096499050505792));
		adl.add( AtomDeriv( AtomID( 9, 16), -0.4165857912449045,0.6196806104791652,-0.7253673232985219,0.1388042638311867,0.06633621771254374,-0.02304572545572723));

		adv.validate_atom_deriv_list( adl );

	}

	/// This can be used to ensure norm matches norm-numeric
	void test_paapp_deriv_check()
	{
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( p_aa_pp, 0.5 );
		MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void dont_test_write_paapp_energies() {
		using namespace core::scoring;
		P_AA const & paa = core::scoring::ScoringManager::get_instance()->get_P_AA();
		// lets output a table of values from the leucine phi/psi distribution
		for ( Size ii = 0; ii < 360; ++ii ) {
			Real ii_phi = -180.0 + ii;
			for ( Size jj = 0; jj < 360; ++jj ) {
				Real jj_psi = -180.0 + jj;
				std::cout << paa.P_AA_pp_energy( core::chemical::aa_leu, ii_phi, jj_psi );
				if ( jj != 359 ) std::cout << " ";
			}
			std::cout << "\n";
		}
	}

};


