// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/membrane_symdocking.cc
///
/// @brief  Membrane Framework Application: Symmetric Protein-Protein Docking in Membranes
/// @details The membrane protien symmetric docking application creates assemblies of c-symmetric
///             complexes, embeds the protein in the membrane, and uses the RosettaMP framework
///             along with RosettaSymDock to model the final complex.
///    Last Modified: 2/9/15
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Project Headers
#include <protocols/symmetric_docking/membrane/MPSymDockMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR( "apps.pilot.membrane.membrane_symdocking" );

/// @brief Main method
int
main( int argc, char * argv [] )
{
	using namespace protocols::jd2;

	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Setup Membrane Symdocking & go!
		using namespace protocols::symmetric_docking::membrane;
		MPSymDockMoverOP mpsymdock( new MPSymDockMover() );
		protocols::jd2::JobDistributor::get_instance()->go( mpsymdock );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

