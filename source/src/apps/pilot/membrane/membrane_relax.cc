// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/membrane_relax.cc
///
/// @brief  Membrane Framework Application: Membrane Protein Structure Refinement
/// @details The membrane framework currently supports membrane protein structure refinement:
///    Combines standard relax cycles with refinement of membrane protein structures
///    using membrane minimization techniques.
///    Last Modified: 7/26/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// App Headers
#include <devel/init.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "apps.pilot.membrane.membrane_relax" );

using namespace protocols::moves;

/// @brief Membrane Relax Main
int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::moves;
		using namespace protocols::membrane;
		using namespace protocols::relax;

		// Initialize Options System, factories, etc
		devel::init( argc, argv );
		protocols::jd2::register_options();

		// Create a new sequence Mover
		SequenceMoverOP seqmov( new SequenceMover() );

		// Create membrane movers & relax mover
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		MembranePositionFromTopologyMoverOP init_pos( new MembranePositionFromTopologyMover() );
		RelaxProtocolBaseOP relax_protocol = generate_relax_from_cmd();

		// Add movers in sequence for memrbane relax protocol
		seqmov->add_mover( add_memb );
		seqmov->add_mover( init_pos );
		seqmov->add_mover( relax_protocol );

		// Execute membrane relax protocol
		protocols::jd2::JobDistributor::get_instance()->go( seqmov );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
