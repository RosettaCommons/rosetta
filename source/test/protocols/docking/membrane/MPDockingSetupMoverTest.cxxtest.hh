// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/docking/membrane/MPDockingSetupMoverTest.cxxtest.hh
/// @brief  Unit test for MPDockingSetupMover, currently only testing parse_my_tag() input
/// @author Aleexsan Adal (aleex.adal@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

static basic::Tracer TR("MPDockingSetupMoverTest");

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

class MPDockingSetupMoverTest : public CxxTest::TestSuite {

public: // test functions

	/// Test Setup Functions ////////

	/// @brief Setup Test
	void setUp(){

		// Initialize
		core_init();

		// Load in pose from pdb
		pose_ = core::pose::PoseOP( new Pose() );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1AFO_AB.pdb" , core::import_pose::PDB_file);

		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/geometry/1AFO_AB.span";

		// AddMembraneMover
		AddMembraneMoverOP add_mem( new AddMembraneMover( spanfile ) );
		add_mem->apply( *pose_ );

	}

	/// @brief Standard Tear Down
	void tearDown() {}


	///// Test Methods /////////////

	////////////////////////////////////////////////////////////////////////////////

	/// @brief test parse_my_tag() read-in
	void test_parse_my_tag_input() {
		std::stringstream tag_ss("<MPDockingSetupMover optimize1=true optimize2=false pose1=\"protocols/membrane/1AFO_AB.pdb\" pose2=\"protocols/membrane/1AFO_AB.pdb\" span1=\"s1.span\" span2=\"s2.span\" />");
		utility::tag::TagCOP tag = utility::tag::Tag::create( tag_ss );

		protocols::docking::membrane::MPDockingSetupMoverOP xmpd(new protocols::docking::membrane::MPDockingSetupMover());
		xmpd->parse_my_tag( tag );

		TS_ASSERT_EQUALS(xmpd->get_optimize1(), true);
		TS_ASSERT_EQUALS(xmpd->get_optimize2(), false);
		TS_ASSERT( (xmpd->get_poses())[1] != NULL);
		TS_ASSERT( (xmpd->get_poses())[2] != NULL);
		TS_ASSERT( !(xmpd->get_spanfiles())[1].compare("s1.span") );
		TS_ASSERT( !(xmpd->get_spanfiles())[2].compare("s2.span") );
	}


private: // data

	core::pose::PoseOP pose_;

};
