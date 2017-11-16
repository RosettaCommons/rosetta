// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Unit headers
#include <devel/domain_assembly/RedesignDomainAssemblyMover.hh>

// Project headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <devel/init.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.andrew.domain_assembly_jd2" );

int main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );

		devel::domain_assembly::DomainAssemblyMoverOP mover( new devel::domain_assembly::DomainAssemblyMover );

		// GO!
		protocols::jd2::JobDistributor::get_instance()->go( mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

