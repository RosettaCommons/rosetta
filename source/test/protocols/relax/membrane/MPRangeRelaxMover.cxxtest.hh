// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/relax/membrane/MPRangeRelaxMover.cxxtest.hh
/// @brief  Unit test for MPRangeRelaxMover, currently only testing parse_my_tag()
/// @author aleex-adal (aleex-adal@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/moves/Mover.hh>

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

static basic::Tracer TR("protocols.relax.membrane.MPRangeRelaxMover.cxxtest");

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

class MPRangeRelaxMover : public CxxTest::TestSuite {

public:

	/// @brief Setup Test
	void setUp(){

		// Initialize
		core_init();

		// Load in pose from pdb
		pose_ = core::pose::PoseOP( new Pose() );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1AFO_AB.pdb", core::import_pose::PDB_file);

		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/geometry/1AFO_AB.span";

		// AddMembraneMover
		AddMembraneMoverOP add_mem( new AddMembraneMover( spanfile ) );
		add_mem->apply( *pose_ );
	}

	/// @brief Standard Tear Down
	void tearDown(){

	}

	///// Test Methods /////////////

	////////////////////////////////////////////////////////////////////////////////

	/// @brief test parse_my_tag() read-in
	void test_mprangerelaxmover_parse_my_tag_input() {
		std::stringstream tag_ss("<MPRangeRelaxMover native=\"protocols/relax/membrane/2mpn_tr_native.pdb\" sfxn=\"mpframework_smooth_fa_2012.wts\" center_resnumber=20 set_tm_helical=true optmem=true/>");
		utility::tag::TagCOP tag = utility::tag::Tag::create( tag_ss );
		protocols::relax::membrane::MPRangeRelaxMoverOP xmprr(new protocols::relax::membrane::MPRangeRelaxMover());

		xmprr->parse_my_tag( tag );

		TS_ASSERT(xmprr->get_native() != nullptr);
		TS_ASSERT(xmprr->get_sfxn() != nullptr);
		TS_ASSERT_EQUALS(xmprr->get_center_resnumber(), 20);
		TS_ASSERT_EQUALS(xmprr->get_set_tm_helical(), true);
		TS_ASSERT_EQUALS(xmprr->get_optmem(), true);
	}

private:
	core::pose::PoseOP pose_;

};
