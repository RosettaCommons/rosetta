// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file 		src/apps/pilot/ralford/transform_into_mem.cc
///
/// @brief 		needs to update docs on this one
/// @details	needs to update docs on this one
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/18/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/TransformIntoMembraneMover.hh> 
#include <protocols/membrane/MembranePositionFromTopologyMover.hh> 
#include <protocols/membrane/visualize/VisualizeMembraneMover.hh> 
#include <protocols/membrane/visualize/ShowMembranePlanesMover.hh> 
#include <protocols/membrane/SetMembranePositionMover.hh> 

#include <core/kinematics/Stub.hh> 
#include <core/kinematics/Jump.hh> 

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh> 
#include <core/kinematics/Edge.hh> 

// Project Headers
#include <protocols/jd2/JobDistributor.hh> 
#include <protocols/jd2/util.hh>

#include <core/pose/Pose.hh> 
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyz.functions.hh> 
#include <numeric/xyzVector.hh> 

#include <numeric/random/random.hh> 
#include <utility/pointer/owning_ptr.hh> 

#include <basic/Tracer.hh> 

// C++ headers
#include <iostream>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.visualize_membrane" );

using namespace protocols::moves;

/// @brief View membrane planes based on normal/center
class MembraneViewMover : public Mover {

public: 

	MembraneViewMover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneViewMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace protocols::membrane::visualize; 
		using namespace numeric; 
		using namespace core::kinematics;

		// Add a Membrane
		AddMembraneMoverOP add_memb = new AddMembraneMover(); 
		add_memb->apply(pose); 

		// Rotate a bit

		//Vector center(16.6831, 3.8529, -1.39162);
		//Vector normal(0.0437651, 0.00068703, 0.999042);

		Vector center( 0, 5, 0 ); //= pose.conformation().membrane_info()->membrane_center(pose.conformation());
		Vector normal( 0, 0, 10 );
		//SetMembranePositionMoverOP transform_memb2 = new SetMembranePositionMover( center, normal, 1 ); 
		//transform_memb2->apply( pose ); 

 
		SetMembranePositionMoverOP transform_memb = new SetMembranePositionMover( center, normal, 2 ); 
		transform_memb->apply( pose ); 


		//MembranePositionFromTopologyMoverOP initialize_memb = new MembranePositionFromTopologyMover(); 
		//initialize_memb->apply( pose ); 

	
	//	Vector normal = pose.conformation().membrane_info()->membrane_normal(pose.conformation());

		//SetMembranePositionMoverOP transform_memb = new SetMembranePositionMover( center, normal, 1 ); 
		//TransformIntoMembraneMoverOP transform_memb = new TransformIntoMembraneMover( 2, center, normal );
		//transform_memb->apply( pose ); 


		pose.conformation().membrane_info()->show();

	}
};

typedef utility::pointer::owning_ptr< MembraneViewMover > MembraneViewMoverOP; 
typedef utility::pointer::owning_ptr< MembraneViewMover const > MembraneViewMoverCOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Minimize with mp
		MembraneViewMoverOP mpview = new MembraneViewMover();
		protocols::jd2::JobDistributor::get_instance()->go( mpview );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

