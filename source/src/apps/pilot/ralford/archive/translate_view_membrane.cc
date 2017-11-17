// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/ralford/translate_view_membrane.cc
///
/// @brief   Translate Structures into a Single Membrane and View
/// @details Given a list of input structures, translate all of those structures
///    into a single membrane. Add a series of membrane virtual residues
///    to construct a set of planes which enable visualization of the membrane
///    posiiton.
///
/// @author  Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/19/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/visualize/VisualizeMembraneMover.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/random/random.hh>
#include <utility/pointer/owning_ptr.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR( "apps.pilot.ralford.translate_view_membrane" );

using namespace protocols::moves;

/// @brief Transform Pose into Membrane and View Membrane
class MembraneViewMover : public Mover {

public:

	MembraneViewMover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneViewMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace protocols::membrane::visualize;

		// Add membrane to pose
		AddMembraneMoverOP add_memb = new AddMembraneMover();
		add_memb->apply( pose );

		// Translate into membrane
		TransformIntoMembraneMoverOP transform_memb = new TransformIntoMembraneMover();
		transform_memb->apply( pose );

		// Add membrane residues to pose to visualize membrane planes
		VisualizeMembraneMoverOP visualize_memb = new VisualizeMembraneMover();
		visualize_memb->apply( pose );

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

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

