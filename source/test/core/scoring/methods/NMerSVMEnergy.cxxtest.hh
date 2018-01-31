// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/NMerSVMEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::NMerSVMEenergy.cc
/// @author Indigo King (indigo.c.king@gmail.com)
/// @author

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/NMerSVMEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>

#include <core/kinematics/DomainMap.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::id;
using namespace core::kinematics;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class NMerSVMEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using namespace core;
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	// --------------- Test Cases --------------- //
	void test_eval_energy()
	{
		Pose trpcage( create_trpcage_ideal_pose() );

		utility::vector1< std::string > svm_fname_vec;
		svm_fname_vec.push_back( "sequence/mhc_svms/HLA-DRB10101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" );
		utility::vector1< std::string > pssm_fname_vec;
		pssm_fname_vec.push_back( "sequence/mhc_pssms/HLA-DRB10101_nooverlap.9mer.norm.pssm" );
		NMerSVMEnergy nmer_svm_energy(
			core::Size( 9 ),
			false,
			core::Size( 3 ),
			true,
			0.0,
			svm_fname_vec,
			pssm_fname_vec
		);

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 1 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], -0.0240031, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 2 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.243551, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 3 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.260367, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 4 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.446798, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 5 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.319999, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 6 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.470622, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 7 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.325501, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 8 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.0427313, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 9 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], -0.0288486, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 10 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], -0.0651805, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 11 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], -0.267079, 1e-4 );
		}

		{
			EnergyMap emap;
			nmer_svm_energy.residue_energy( trpcage.residue( 12 ), trpcage, emap );
			TS_ASSERT_DELTA( emap[ nmer_svm ], 0.165212, 1e-4 );
		}

	}


};


