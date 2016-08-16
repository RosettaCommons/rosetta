// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/mp_viewer.cc
///
/// @brief  RosettaMP Real time visualization of membrane proteins
/// @details Use the PyMOL viewer to visualize the membrane planes based on position
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @note   Last Updated: 5/18/15

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/moves/PyMolMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

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

/// @brief Quick Container Mover: Visualize Membrane Protein Using the PyMol Viewer
class ViewMembraneProteinMover : public Mover {

public:

	/// @brief Default Constructor
	ViewMembraneProteinMover() :
		Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "ViewMembraneProteinMover"; }

	/// @brief Setup the Membrane Framework. Flags should take care of the rest
	/// of visualization
	void apply( Pose & pose ) {

		using namespace basic::options;
		using namespace protocols::membrane;
		using namespace protocols::moves;

		// Add Membrane
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		add_memb->apply( pose );

		// Send a set of viewable planes to pymol
		PyMolMoverOP pymol_mover( new PyMolMover() );
		pymol_mover->apply( pose );

		if ( option[ OptionKeys::mp::setup::position_from_topo ].user() ) {

			// Initialize Membrane
			MembranePositionFromTopologyMoverOP init_memb( new MembranePositionFromTopologyMover() );
			init_memb->apply( pose );
		}
	}
};

// Hook Mover into Rosetta
typedef utility::pointer::shared_ptr< ViewMembraneProteinMover > ViewMembraneProteinMoverOP;
typedef utility::pointer::shared_ptr< ViewMembraneProteinMover const > ViewMembraneProteinMoverCOP;

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

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
