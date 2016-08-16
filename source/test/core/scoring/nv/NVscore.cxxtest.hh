// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/nv/NVScore.cxxtest.hh
/// @brief  test suite for core::scoring::NVScore.cc
/// @author Sam DeLuca

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/nv/NVscore.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>


//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::nv;

class NVscoreTests : public CxxTest::TestSuite {

public:

	PoseOP the_pose;
	NVscoreOP nv_score;

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

		nv_score = NVscoreOP( new NVscore );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
		nv_score.reset();
	}


	// --------------- Test Cases --------------- //
	void test_eval_energy() {

		float const TOLERATED_ERROR = 0.0001;

		Real neigh_vect_answers[] = {
			-0.648969, 0.391302, 2.69084, -0.359112, -1.01807, -0.213835, -0.25529, -1.20491, -1.39418, -0.716153,
			-0.352967, -1.76049, -0.319042, -0.344338, -0.380177, -1.41007, -0.896221, -0.361165, -0.288861, 0.287146,
			0.129568, 0.401269, -0.348142, -0.510714, 0.135269, -0.992935, -0.115201, 0.0922533, -0.193909, -0.127102,
			0.592978, -0.0962339, -0.60012, 1.8265, -0.31043, -1.54499, -0.424731, -0.63256, -0.66194, -1.70508, -0.706986,
			-0.316944, -0.27182, -0.340076, 0.0030844, -0.531556, -0.378463, -0.33029, -0.0955399, -0.527367, -1.29255, -1.08097,
			-0.724529, -0.359241, -0.44333, 0.542265, -0.509971, -0.926208, -0.283379, -0.933063, -1.49854, -0.5465, -0.202579,
			-1.14428, -0.678482, -0.516606, -0.294236, 0.385786, -0.311551, -0.303241, -0.773431, 0.68871, 0.487456, -1.67176,
			-1.1281, -0.749384, -0.7242, -0.218217, -0.317773, -0.758513, 0.173812, -0.393873, -0.627719, -0.0492903, -0.359348,
			-1.63023, -0.358514, -0.553371, -0.349237, -0.0244665, -0.669393, -0.181687, -0.678182, -1.13211, 0.235373, -0.00100187,
			-0.276021, 0.454144, 0.0296966, -0.858374, -1.52818, -0.671063, 1.41223, -0.362511, 0.748328, 1.88056, -0.448444, -0.17265,
			-0.470131, -0.384258, -0.679863, -0.629524, 0.021603};


		Real neigh_vect_raw_answers [] = {
			0.476311, 0.66641, 0.707215, 0.438855, 0.0785003, 0.472318, 0.437379, 0.0841765, 0.153148, 0.507288, 0.340164, 0.10929,
			0.340478, 0.605007, 0.290884, 0.0964949, 0.510911, 0.552303, 0.34806, 0.421362, 0.0578164, 0.453584, 0.363073, 0.582907,
			0.37978, 0.10756, 0.247582, 0.605969, 0.379011, 0.459564, 0.49967, 0.611817, 0.347269, 0.660863, 0.427586, 0.102155,
			0.289371, 0.511263, 0.246698, 0.112988, 0.415312, 0.641825, 0.299329, 0.35094, 0.658633, 0.342121, 0.501493, 0.0632692,
			0.334978, 0.554896, 0.170707, 0.0804143, 0.355673, 0.248904, 0.154644, 0.112302, 0.0719913, 0.137312, 0.305023, 0.0311858,
			0.11239, 0.404279, 0.23968, 0.088159, 0.401521, 0.651047, 0.457991, 0.724598, 0.330538, 0.573777, 0.215823, 0.516548, 0.473279,
			0.126669, 0.0866042, 0.552568, 0.370117, 0.0332649, 0.240645, 0.243362, 0.721868, 0.539603, 0.434351, 0.640276, 0.441751, 0.107502,
			0.314473, 0.526529, 0.325808, 0.0551826, 0.396212, 0.61032, 0.286425, 0.194219, 0.72817, 0.349774, 0.412975, 0.068646, 0.541069, 0.522818,
			0.12849, 0.32179, 0.663596, 0.31052, 0.764365, 0.723574, 0.147289, 0.296004, 0.602453, 0.505124, 0.438797, 0.261803, 0.372858,
			};

		EnergyMap emap;
		the_pose->update_residue_neighbors();
		for ( int ii = 1; ii <= 113; ++ii ) {
			emap.zero();
			nv_score->residue_energy( the_pose->residue(ii), *the_pose, emap );


			TS_ASSERT_DELTA( emap[ neigh_vect ], neigh_vect_answers[ ii - 1 ], TOLERATED_ERROR );
			TS_ASSERT_DELTA( emap[ neigh_vect_raw], neigh_vect_raw_answers[ii - 1], TOLERATED_ERROR);
		}

	}

};


