// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/FaMbEnv.cxxtest.hh
///
/// @brief   Unit Test: FaMbEnv Object
/// @details Unit test of the FaMbEnv score term
///         from the Rosetta Membrane score function.
///      Includes a test to make sure the analytical
///      and numeric derivatives match.
///    Last Modified: 5/27/2021
///
/// @author  Hope Woods (hope.woods@vanderbilt.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>


//Unit Headers

// Package Headers


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>


#include <basic/Tracer.hh>


// Utility Headers

using namespace core;
using namespace core::conformation;

static basic::Tracer TR("core.scoring.membrane.FaMbEnv.cxxtest");

class FaMbEnvTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {
		using namespace core::import_pose;
		using namespace core::pose;


		// Initialize core & options system
		core_init_with_additional_options( "-in:file:spanfile core/conformation/membrane/1AFO_AB.span" );
		TR << "core_init" << std::endl;

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods /////////////////////////////////

	void test_fambenvenergy_deriv() {
#ifndef MULTI_THREADED //the membrane score terms are non-threadsafe
		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;


		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_mbenv, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6);
#endif

	}

private:



}; // FaMbEnv unit test
