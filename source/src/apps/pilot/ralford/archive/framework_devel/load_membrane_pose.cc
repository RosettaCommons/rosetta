// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file 		src/apps/pilot/ralford/load_membrane_pose.cc
///
/// @brief 		Rosetta Membrane Framework - Bottom Level Integration Test
/// @details    Check that a membrane pose is being loaded in reliably, coordinates make
///				sense, and the basic object machinery (SpanningTopology, MembraneInfo, etc)
///				are working just fine
///
/// 			Last Modified: 6/28/14
///				Note: Part of the Membrane Framework Applications Suite
///				Versioning: 1.0
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/MembranePositionFromTopologyMover.hh> 

#include <protocols/jd2/JobDistributor.hh> 
#include <protocols/jd2/util.hh>

#include <protocols/moves/Mover.hh> 

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 

#include <core/pose/Pose.hh> 
#include <core/types.hh> 

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>
#include <cstdlib> 

using namespace protocols::moves; 

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using Rosetta's new Membrane Framework
class LoadMembraneMover : public Mover {

public: 

	/// @brief Default Constructor
	LoadMembraneMover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "LoadMembraneMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace protocols::moves; 

		// Add Membrane
		AddMembraneMoverOP add_memb( new AddMembraneMover() ); 
		add_memb->apply( pose ); 

		// Show new pose membrane
		pose.conformation().membrane_info()->show();
		
		// Initialize Membrane
		MembranePositionFromTopologyMoverOP init_memb( new MembranePositionFromTopologyMover() );
		init_memb->apply( pose ); 

		// Show new pose membrane
		pose.conformation().membrane_info()->show();

	}
};

typedef utility::pointer::shared_ptr< LoadMembraneMover > LoadMembraneMoverOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		using namespace protocols::moves;
		using namespace protocols::membrane;

		protocols::jd2::register_options();

		// Create and kick off a new load membrane mover
		LoadMembraneMoverOP load_memb( new LoadMembraneMover() ); 
		protocols::jd2::JobDistributor::get_instance()->go( load_memb );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

