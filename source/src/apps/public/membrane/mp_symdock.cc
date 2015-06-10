// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		apps/pilot/membrane/mp_symdock.cc
///
/// @brief		RosettaMP Membrane Symmetric Protien-Protein Docking Application
/// @details	Assemble C-Symmetric poses in the membrane
///
/// @author 	Rebecca Faye Alford (rfalford12@gmail.com)
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @note 		Last Updated: 5/18/15

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

static basic::Tracer TR( "apps.public.membrane.mp_symdock" );

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

