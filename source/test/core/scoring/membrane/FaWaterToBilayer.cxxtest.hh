// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/FaWaterToBilayer.cxxtest.hh
///
/// @brief   Unit Test: FaWaterToBilayer Object
/// @details This suite is for the FaWaterToBilayer energy term.
///          It includes a test to make sure the analytical
///          and numeric derivatives are equal.
///
/// @author  Hope Woods (hope.woods@vanderbilt.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>


//Unit Headers
#include <protocols/membrane/scoring/FaWaterToBilayerEnergy.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

// Utility Headers
#include <utility/vector1.hh>

using namespace core;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;

static basic::Tracer TR("core.scoring.membrane.FaWaterToBilayer.cxxtest");

class FaWaterToBilayerTest : public CxxTest::TestSuite {

public: // test functions

	// Test Setup Functions ///////////////////////////

	/// @brief Setup Test
	void setUp() {
		using namespace core::import_pose;
		using namespace core::pose;



		core_init_with_additional_options( "-in:membrane" );
		// Initialize core & options system
		TR << "core_init" << std::endl;

	}

	/// @brief Tear Down Test
	void tearDown() {}

	// Test Methods ////////////////////////////////


	void test_fampenvenergy_deriv_1py6_pore() {
		using namespace core::scoring;
		using namespace core::import_pose;
		using namespace core::pose;

		// Load in pose from pdb
		core::pose::PoseOP pose = utility::pointer::make_shared< Pose >();
		TR << "initialize pose" << std::endl;
		core::import_pose::pose_from_file( *pose, "core/conformation/membrane/1PY6_mp_coords.pdb", core::import_pose::PDB_file);
		TR << "pose from file" << std::endl;

		// Initialize Spans from spanfile
		std::string span = "core/conformation/membrane/1PY6.span";

		AddMembraneMoverOP add_memb( new AddMembraneMover( span ));

		add_memb->apply( *pose );
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_water_to_bilayer, 0.3 );
		core::kinematics::MoveMap movemap;
		movemap.set_bb( true );
		AtomDerivValidator adv( *pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5);

	}



private:



}; // FaWaterToBilayer unit test
