// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/OptimizeProteinEmbeddingMoverTest.cxxtest.hh
/// @brief  Unit test for OptimizeProteinEmbeddingMover, currently only testing parse_my_tag()
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
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
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
#include <basic/datacache/DataMap.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR("OptimizeProteinEmbeddingMoverTest");

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

class OptimizeProteinEmbeddingMoverTest : public CxxTest::TestSuite {

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


	/// @brief test parse_my_tag() read-in
	//TODO: modify this for future options that may be added
	void test_parse_my_tag_input() {
		std::stringstream tag_ss("<OptimizeProteinEmbeddingMover option=option_value/>");
		utility::tag::TagCOP tag = utility::tag::Tag::create( tag_ss );

		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		protocols::moves::Movers_map mm;
		core::pose::Pose pose;

		OptimizeProteinEmbeddingMoverOP xpro( new OptimizeProteinEmbeddingMover() );
		xpro->parse_my_tag( tag, dm, fm, mm, pose );

		TS_ASSERT(true);
	}


private: // data
	core::pose::PoseOP pose_;

};
