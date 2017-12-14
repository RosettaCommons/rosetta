// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/view_membrane_protein.cc
///
/// @brief  Membrane Framework Application: Real-Time Visualization of Membrane Proteins
/// @details The PyMOL mover currently includes methods for visualizing the membrane
///    planes. This can be turned on using the flags -show_simulation_in_pymol
///    Last Modified: 7/26/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/moves/PyMOLMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// C++ headers
#include <iostream>
#include <cstdlib>

using namespace protocols::moves;

/// @brief Quick Container Mover: Visualize Membrane Protein Using the PyMOL Viewer
class ViewMembraneProteinMover : public Mover {

public:

	/// @brief Default Constructor
	ViewMembraneProteinMover() :
		Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const override { return "ViewMembraneProteinMover"; }

	/// @brief Setup the Membrane Framework. Flags should take care of the rest
	/// of visualization
	void apply( Pose & pose ) override {

		using namespace basic::options;
		using namespace protocols::membrane;
		using namespace protocols::moves;

		// Add Membrane
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		add_memb->apply( pose );

		// Send a set of viewable planes to pymol
		PyMOLMoverOP pymol_mover( new PyMOLMover() );
		pymol_mover->apply( pose );

		if ( option[ OptionKeys::mp::setup::position_from_topo ].user() ) {

			// Initialize Membrane
			MembranePositionFromTopologyMoverOP init_memb( new MembranePositionFromTopologyMover() );
			init_memb->apply( pose );
		}
	}
};

// Hook Mover into Rosetta
using ViewMembraneProteinMoverOP = utility::pointer::shared_ptr<ViewMembraneProteinMover>;
using ViewMembraneProteinMoverCOP = utility::pointer::shared_ptr<const ViewMembraneProteinMover>;

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
		ViewMembraneProteinMoverOP view_memb( new ViewMembraneProteinMover() );
		protocols::jd2::JobDistributor::get_instance()->go( view_memb );

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
