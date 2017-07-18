// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/ralford/membrane_relax.cc
///
/// @brief   Membrane Relax Application
/// @details    High resolution FastRelax using a custom set of settings & modifications
///    for refinement of membrane protein structure. Adaptation of the current FastRelax
///    protocol
///    Last Modified: 7/21/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>

// Project Headers
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/relax_main.hh>
#include <protocols/relax/util.hh>

// Project Headers
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

using namespace protocols::moves;
using namespace protocols::jd2;

/// @brief Membrane Relax Mover: Create a Membrane Pose and Apply Relax Protocol
class MembraneRelaxMover : public Mover {

public:

	MembraneRelaxMover() :
		Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneRelaxMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace protocols::membrane;
		using namespace protocols::moves;
		using namespace protocols::relax;

		// Create a membrane pose
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		add_memb->apply( pose );

		// Initialze the posiiton of a membrane protein
		MembranePositionFromTopologyMoverOP init_position( new MembranePositionFromTopologyMover() );
		init_position->apply( pose );

		// Setup Relax Base Protocol
		protocols::moves::MoverOP protocol = generate_relax_from_cmd();
		protocols::jd2::set_native_in_mover( *protocol );
		protocol->apply(pose);

	}
};

typedef utility::pointer::shared_ptr< MembraneRelaxMover > MembraneRelaxMoverOP;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::moves;
		using namespace scoring;
		using namespace basic::options;
		using namespace protocols::membrane;

		protocols::relax::ClassicRelax::register_options();
		protocols::jd2::register_options();

		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::in::file::movemap );
		option.add_relevant( OptionKeys::relax::fast );

		devel::init(argc, argv);

		// Create and kick off a new relax mover
		MembraneRelaxMoverOP mprlx( new MembraneRelaxMover() );
		protocols::jd2::JobDistributor::get_instance()->go( mprlx );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

