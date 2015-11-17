// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/P_AA_pp_Energy_symmgly.cxxtest.hh
/// @brief  test suite for the P_AA_pp score term with the -symmetric_gly_tables option turned on.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <basic/Tracer.hh>
#include <core/scoring/methods/P_AA_pp_Energy.hh>
#include <core/scoring/P_AA.hh>
#include <core/scoring/ScoringManager.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>


#include <numeric/conversions.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

static basic::Tracer TR("core.scoring.methods.P_AA_pp_Energy_symgly.cxxtest");

class P_AA_pp_Energy_SymGlyTests : public CxxTest::TestSuite {

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
		core_init_with_additional_options("-score:symmetric_gly_tables");

		paapp_energy = P_AA_pp_EnergyOP( new P_AA_pp_Energy );
		the_pose = core::pose::PoseOP(new core::pose::Pose);

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
		paapp_energy.reset();
	}


	// --------------- Test Cases --------------- //
	void test_eval_symm_energy()
	{
		core::pose::make_pose_from_sequence(*the_pose, "GGGGGG", "fa_standard", false);
		the_pose->set_phi(2, -61);
		the_pose->set_psi(2, -41);
		the_pose->set_phi(3, 61);
		the_pose->set_psi(3, 41);
		the_pose->set_phi(4, -45);
		the_pose->set_psi(4, 120);
		the_pose->set_phi(5, 45);
		the_pose->set_psi(5, -120);

		float const TOLERATED_ERROR = 0.0001;

		EnergyMap emap2, emap3, emap4, emap5;

		paapp_energy->residue_energy( the_pose->residue(2), *the_pose, emap2 );
		paapp_energy->residue_energy( the_pose->residue(3), *the_pose, emap3 );
		paapp_energy->residue_energy( the_pose->residue(4), *the_pose, emap4 );
		paapp_energy->residue_energy( the_pose->residue(5), *the_pose, emap5 );

		TR << "Energy: " << emap2[p_aa_pp] << " Mirror: " << emap3[p_aa_pp] << std::endl;
		TR << "Energy: " << emap4[p_aa_pp] << " Mirror: " << emap5[p_aa_pp] << std::endl;
		TR.flush();

		TS_ASSERT_DELTA( emap2[p_aa_pp], emap3[p_aa_pp], TOLERATED_ERROR );
		TS_ASSERT_DELTA( emap4[p_aa_pp], emap5[p_aa_pp], TOLERATED_ERROR );

		return;
	}

};


