// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    score_membrane.cc
///
/// @brief   Integration Test Precursor - Score a Membrane Pose
/// @details Score a membrane pose within the membrane protein framework
///			 Last Modified: 3/11/14
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/CreateMembranePoseMover.hh>

#include <core/conformation/membrane/SpanningTopology.hh> 
#include <core/membrane/geometry/EmbeddingFactory.hh> 

// Package Headers
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core::conformation::membrane;
using namespace core::membrane::geometry;

static basic::Tracer TR( "apps.pilot.ralford.embed_mp_pose_from_topo" );

void*
my_main( void* )
{
	using namespace protocols::membrane;
	using namespace protocols::jd2;
	
	// Create a Membrane Pose
	CreateMembranePoseMoverOP mp = new CreateMembranePoseMover();
	JobDistributor::get_instance()->go(mp);
	core::pose::PoseOP pose = mp->get_membrane_pose();
	
	// Get Mmebrane Topology for the Pose and create a dummy definiiton
	utility::vector1< SpanningTopology > topology = pose->conformation().membrane_info()->spanning_topology();
	EmbeddingFactoryOP emb = new EmbeddingFactory();

	// Compute Embedding for Chain A
	EmbedConfigInfoOP embA = emb->embed_from_topology( pose, topology[1] );

	// Compute Embedding for Chain B
	EmbedConfigInfoOP embB = emb->embed_from_topology( pose, topology[2] );

	// Print Resulting Embeddings for Chain B
	TR << "Calculated Membrane Topology Based Embedding for Pose Chain A" << std::endl;
	TR << "Normal: " << embA->normal.x() << " " << embA->normal.y() << " " << embA->normal.z() << " " << std::endl;
	TR << "Center: " << embA->center.x() << " " << embA->center.y() << " " << embA->center.z() << " " << std::endl;
	
	return NULL;
}

///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
}

