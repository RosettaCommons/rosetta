// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MMTorsionScore.cxxtest.hh
/// @brief  test suite for core:scoring::methods::MMTorsionEnergy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/MMTorsionEnergy.hh>
#include <core/types.hh>

// Project headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
// Auto-header: duplicate removed #include <core/io/pdb/pose_io.hh>

// Utility headers

// C++ headers
#include <iostream>
// AUTO-REMOVED #include <iomanip>

//Auto Headers
#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

// --------------- Test Class --------------- //

class MMTorsionEnergyTests : public CxxTest::TestSuite {

public:

	PoseOP pose;
	MMTorsionEnergyOP mmtorsionenergy;
	Real delta;

	// --------------- Suite-level Fixture --------------- //

	MMTorsionEnergyTests() {
		core_init_with_additional_options( "-no_optH" );
	}

	virtual ~MMTorsionEnergyTests() {}

	static MMTorsionEnergyTests* createSuite() {
		return new MMTorsionEnergyTests();
	}

	static void destroySuite( MMTorsionEnergyTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	void setUp() {
		// init pose
		pose = create_test_in_pdb_poseop();
		//core::import_pose::pose_from_pdb( *pose, "core/scoring/methods/test_in.pdb" );

		//init mmtorsionenergy
		mmtorsionenergy = new MMTorsionEnergy;

		// init delta
		delta = 0.0001;

	}

	void tearDown() {
		pose = 0;
		mmtorsionenergy = 0;
	}

	// --------------- Test Cases --------------- //

	void test_residue_pair_energy() {

		// pre-computed score
		Real scores[] = { 0,
			1.2042,  1.1870,  1.4731,  1.2514,  1.1907,  1.6506,  1.1967,  1.5185,  1.0828,  1.1470,  1.2589,  1.1156,  1.4611,
			3.3257,  1.8222,  1.0697,  4.3641,  3.3269,  1.8179,  2.4757,  1.2369,  0.2433,  1.2291,  0.9515,  1.3583,  1.2904,
			1.2506,  1.3435,  0.8344,  0.9166,  0.3086,  0.4203,  1.1993,  1.2781,  1.3260,  1.1519,  1.1339,  1.0947,  1.0376,
			1.2534,  1.1034,  1.4087,  2.3028,  1.4376,  2.3747,  0.4104,  0.2811,  1.2602,  1.2723,  1.4752,  1.1245,  1.1916,
			1.1928,  1.3221,  1.2586,  1.2636,  1.2450,  1.2768,  1.3310,  1.3022,  1.2627,  1.4303,  1.2576,  2.1275,  1.9925,
			3.8837,  1.2714,  1.3105,  2.8145,  1.6161,  1.5188,  1.2615,  1.4074,  1.2392,  1.0521,  1.2975,  2.0275,  2.5030,
			1.7383,  0.3880,  4.3318,  0.1847,  1.1258,  1.1631,  1.2481,  1.4338,  1.1627,  1.2471,  1.2130,  1.1967,  1.2531,
			1.7542,  2.2779,  4.2229,  3.2654,  2.5855,  3.4270,  3.3056,  1.1800,  1.2188,  1.2549,  1.8242,  2.8540,  0.6088,
			1.1213,  1.3004,  0.6941,  0.2171,  1.3558,  3.4008,  1.3599,  1.4830,  1.2317,  7.0035,  3.8738 };

		ScoreFunction sfxn; // unused
		EnergyMap emap;
		for ( int i = 1; i <=116; ++i )
			{
			if((i+1) <= 116){
				emap.zero();
				mmtorsionenergy->residue_pair_energy( pose->residue(i), pose->residue(i+1), *pose, sfxn, emap);
				//std::cout <<   std::setprecision(4) << std::fixed << std::setw(8) << emap[ mm_twist ] << ",";
				TS_ASSERT_DELTA( emap[ mm_twist ], scores[i], delta );
			}
			}
	}

	void test_eval_intrares_energy() {

		// pre-computed score
		Real scores[] = {0,
			4.5643,  1.0334,  3.7657,  3.9444,  2.4887,  4.7526,  2.6717,  3.6276,  7.1570,  1.5570,  2.9471,  1.5287,  5.1856,
			5.2445,  2.3671,  6.9671,  7.4941,  5.8335,  6.5351,  5.4304,  6.8873,  3.2784,  4.6727, 10.8303,  2.0761,  3.3922,
			7.4531,  4.6213,  3.3336,  0.0000,  1.7819,  3.7750,  2.5698,  3.0714,  1.1695,  1.4444,  2.0391,  2.4508,  2.0292,
			2.4272,  3.5271, 10.6529,  5.1524,  4.7644,  0.0000,  6.9020,  5.0537,  6.8592,  0.0000,  5.6415,  2.2000,  1.4329,
			7.5779,  1.9759,  3.9103,  1.9197,  3.4024,  4.9650,  5.4159,  1.8099,  1.3007,  4.3203,  6.9618,  2.0340,  9.4489,
			4.7300,  4.7039,  5.0189,  3.7708,  8.4631,  4.0119,  5.4095,  4.8868,  1.8824,  0.8229,  1.2398,  6.8120,  4.3998,
			0.0000,  3.1745,  3.9023,  2.4695,  1.5024,  1.8981,  5.2167,  2.2636,  5.2081,  3.7651,  5.7126,  3.2262,  7.5930,
			1.6236,  2.4714,  3.3233,  2.7632,  2.1545,  5.0748,  7.2648,  3.2494,  3.7799,  2.3941,  7.3457,  6.0559,  6.2063,
			1.6434,  5.9272,  6.8520,  0.0000,  6.4328,  2.7011,  4.6752,  2.8682,  1.9902,  2.3691,  5.6557,  2.9680 };

		ScoreFunction sfxn; // unused
		EnergyMap emap;
		for ( int i = 1; i <=116; ++i )
			{
			emap.zero();
			mmtorsionenergy->eval_intrares_energy( pose->residue(i), *pose, sfxn, emap);
			//std::cout <<   std::setprecision(4) << std::fixed << std::setw(8) << emap[ mm_twist ] << ",";
			TS_ASSERT_DELTA( emap[ mm_twist ], scores[i], delta );
			}
	}
};
