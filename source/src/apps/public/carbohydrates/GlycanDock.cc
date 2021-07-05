// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/carbohydrates/GlycanDock.cc
/// @brief Perform local docking and refinement of a glycoligand into the putative pocket of a protein receptor
/// @author Morgan Nance (morganlnance@gmail.com) based off of dock_glycans app by Jason Labonte (JWLabonte@jhu.edu)


// Unit Headers
#include <devel/init.hh>

// Project Headers
#include <protocols/glycan_docking/GlycanDockProtocol.hh>

#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>


static basic::Tracer TR( "apps.public.carbohydrates.GlycanDock" );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace protocols;

	try {
		devel::init(argc, argv);

		//  protocols::moves::MoverOP glycan_dock_mover
		protocols::glycan_docking::GlycanDockProtocolOP glycan_dock
			( utility::pointer::make_shared
			< protocols::glycan_docking::GlycanDockProtocol >() );

		protocols::jd2::JobDistributor::get_instance()->go(glycan_dock);

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
