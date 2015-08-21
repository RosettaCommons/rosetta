// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/UnfoldedStateEnergy.cxxtest.hh
/// @brief  test suite for UnfoldedStateEnergy energy method
/// @author Ron Jacak

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/UnfoldedStateEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreType.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class UnfoldedStateEnergyTests : public CxxTest::TestSuite {

public:

	PoseOP the_pose;
	Real TOLERATED_ERROR;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		TOLERATED_ERROR = 0.001;
		core_init();

		the_pose = PoseOP( new Pose );
		core::chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::pose::make_pose_from_sequence( *the_pose, "DFGLK", *rsd_set );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
	}


	// --------------- Test Cases --------------- //

	// Test the three different ways of using the UnfoldedStateEnergy energy method class. The difference lies in how the class
	// is constructed: with the default constructor, with an empty EMap, or with a full EMap.

	void test_residue_energy_default_file_weights() {
		// AA fa_atr fa_rep fa_sol fa_intra_rep pro_close fa_pair hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc rama omega fa_dun p_aa_pp
		// WEIGHT 0.8 0.44 0.65 0.004 1.0 0.49  0.585 1.17 1.17 1.1 0.2 0.5 0.56 0.32
		// PHE 3.6792 -0.7291 -1.7340 -6.6803 -0.0053 0.0000 0.4257 0.0216 0.0055 0.0000 0.1227 -0.2181 -1.0325 0.1014
		// GLY 1.6906 -0.4879 -1.2009 -0.0181 -0.0033 0.0000 0.3093 0.0139 0.0093 0.0000 -0.6594 -0.1474 0.0000 1.0637

		UnfoldedStateEnergyOP unweighted_unf_energy( new UnfoldedStateEnergy( scoring::UNFOLDED_SCORE12 ) );

		EnergyMap emap;
		emap.zero();

		unweighted_unf_energy->residue_energy( the_pose->residue(2), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], 1.114, TOLERATED_ERROR );

		emap.zero();
		unweighted_unf_energy->residue_energy( the_pose->residue(3), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], 0.6967, TOLERATED_ERROR );

		unweighted_unf_energy = 0;

	}

	void test_residue_energy_no_weights() {

		EnergyMap emap; // same emap used for initalization and holding scores
		emap.zero();

		UnfoldedStateEnergyOP unweighted_unf_energy( new UnfoldedStateEnergy( scoring::UNFOLDED_SCORE12, emap ) );

		unweighted_unf_energy->residue_energy( the_pose->residue(2), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], 0.0, TOLERATED_ERROR );

		emap.zero();
		unweighted_unf_energy->residue_energy( the_pose->residue(3), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], 0.0, TOLERATED_ERROR );

		unweighted_unf_energy = 0;
	}

	void test_residue_energy_emap_weights() {

		EnergyMap emap;
		emap[ fa_atr ] = 0.8;
		emap[ fa_rep ] = 0.454238;
		emap[ fa_sol ] = 1.04457;
		emap[ fa_intra_rep ] = 0.223082;
		emap[ pro_close ] = 0.374325;
		emap[ fa_pair ] = 0.743825;
		emap[ hbond_sr_bb ] = 0.585;
		emap[ hbond_lr_bb ] = 1.17;
		emap[ hbond_bb_sc ] = 2.58097;
		emap[ hbond_sc ] = 0.477475;
		emap[ dslf_ss_dst ] = 1;
		emap[ dslf_cs_ang ] = 1;
		emap[ dslf_ss_dih ] = 1;
		emap[ dslf_ca_dih ] = 1;
		emap[ rama ] = 0.2;
		emap[ omega ] = 0.151806;
		emap[ fa_dun ] = 0.358969;
		emap[ p_aa_pp ] = 1.01463;
		emap[ unfolded ] = 0.61846;
		emap[ ref ] = 0.0;

		UnfoldedStateEnergyOP weighted_unf_energy_emap( new UnfoldedStateEnergy( scoring::UNFOLDED_SCORE12, emap ) );

		// atr  rep  sol  intra_rep  pro  pair  hb_sr_bb  hb_lr_bb  hb_bb_sc  hb_sc  rama  omega  dun  paapp
		// LEU -3.1780 0.71276 1.75359 2.00079 0.01865 0.00000 0.00000 0.00000 -0.0044 0.00000 -0.1056 0.20147 0.60841 -0.1721
		// LYS 2.8124 -0.6590 -1.8136 -1.0798 -0.0047 0.1160 0.4479 0.0148 0.0091 0.0012 0.2868 -0.1587 -1.5303 0.0155
		// ASP 3.0032 -0.6757 -2.4590 -1.4995 -0.0050 0.1044 0.4009 0.0166 0.0416 0.0125 0.0134 -0.1618 -0.8733 0.1915
		// values below verified using Excel spreadsheet for correctness
		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(4), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], 0.30581, TOLERATED_ERROR );

		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(5), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], -0.29715, TOLERATED_ERROR );

		emap.zero();
		weighted_unf_energy_emap->residue_energy( the_pose->residue(1), *the_pose, emap );
		TS_ASSERT_DELTA( emap[ unfolded ], -0.50548, TOLERATED_ERROR );

		weighted_unf_energy_emap = 0;
	}

};


