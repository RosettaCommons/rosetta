// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/BetaNov16WithAutoSetupMetalsTests.cxxtest.hh
/// @brief  Unit test of -beta_nov16 scorefunction with -auto_setup_metals.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ResidueType.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("BetaNov16WithAutoSetupMetalsTests");


class BetaNov16WithAutoSetupMetalsTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-beta -auto_setup_metals" );

	}

	void tearDown(){

	}


	void test_auto_setup_metals_with_ser_mn(){
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/scoring/4mmx_context_0001_renumbered_trimmed.pdb" , core::import_pose::PDB_file);
		TS_ASSERT_EQUALS( pose.residue_type(12).name3(), "SER" );
		TS_ASSERT( !(pose.residue_type(12).has("HG")) ); //The HG atom should have been stripped by -auto_setup_metals.

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		(*sfxn)(pose);

		core::scoring::ScoreFunctionOP sfxn2( new core::scoring::ScoreFunction );
		sfxn2->set_weight( core::scoring::hxl_tors, 1.0 );
		(*sfxn2)(pose);
	}



};
