// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief   Unit test for SpinAroundPartnerMover
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/SpinAroundPartnerMover.hh>
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

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

class SpinAroundPartnerMoverTest : public CxxTest::TestSuite {

public: // test functions

	/// Test Setup Functions ////////

	/// @brief Setup Test
	void setUp(){

		// Initialize
		core_init();

		// Load in pose from pdb
		pose_ = Pose();
		core::import_pose::pose_from_pdb( pose_, "protocols/membrane/1C17_tr.pdb" );

		// Initialize Spans from spanfile
		std::string spanfile = "protocols/membrane/1C17_tr.span";

		// AddMembraneMover
		AddMembraneMoverOP add_mem( new AddMembraneMover( spanfile ) );
		add_mem->apply( pose_ );

	}

	/// @brief Standard Tear Down
	void tearDown() {}

	///// Test Methods /////////////

	////////////////////////////////////////////////////////////////////////////////

	// test constructor from jumpnum and angle
	void test_constructor_from_jumpnum () {

		TS_TRACE("\n\n========== TESTING CONSTUCTOR FROM JUMP NUMBER");

		// jump = 1
		// angle = 43
		// axis = perpendicular to axis connecting protein embedding centers

		// compute downstream empedding
		SpanningTopologyOP topo = pose_.conformation().membrane_info()->spanning_topology();
		SpanningTopologyOP topo_up_( new SpanningTopology() );
		SpanningTopologyOP topo_down_( new SpanningTopology() );

		// split_topology_by_jump_noshift
		split_topology_by_jump_noshift( pose_, 1, topo, topo_up_, topo_down_ );

		// compute embedding for downstream partner
		EmbeddingDefOP emb_down( compute_structure_based_embedding( pose_, *topo_down_ ) );

		// define vectors
		Vector cntr_before(18.6088, 11.9908, -0.248);
		Vector norm_before(0.0076786, 0.0858815, 0.996276);
		Vector cntr_after(10.0, 12.0, 0.0);
		Vector norm_after(0.0076786, 0.0858815, 0.996276);

		// compare before
		TS_ASSERT( position_equal_within_delta( emb_down->center(), cntr_before, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( emb_down->normal(), norm_before, 0.001 ) );

		// create and run SpinAroundPartnerMover
		SpinAroundPartnerMoverOP spin( new SpinAroundPartnerMover( 1 ) );
		spin->set_x( 10.0 );
		spin->set_y( 12.0 );
		spin->apply( pose_ );

		// compute embedding for downstream partner
		EmbeddingDefOP emb_down1( compute_structure_based_embedding( pose_, *topo_down_ ) );

		// compare after
		TS_ASSERT( position_equal_within_delta( emb_down1->center(), cntr_after, 0.001 ) );
		TS_ASSERT( position_equal_within_delta( emb_down1->normal(), norm_after, 0.001 ) );

	}

	////////////////////////////////////////////////////////////////////////////////

	/// @brief Position equal within delta (helper method)
	bool position_equal_within_delta( Vector a, Vector b, Real delta ) {

		TS_ASSERT_DELTA( a.x(), b.x(), delta );
		TS_ASSERT_DELTA( a.y(), b.y(), delta );
		TS_ASSERT_DELTA( a.z(), b.z(), delta );

		return true;
	}

private: // data

	core::pose::Pose pose_;

};

