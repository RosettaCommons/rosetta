// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Mover partners apart and relax them separately
/// @details Run quick relax on separated partners; this emulates unbound docking
/// @author     JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/QuickRelaxPartnersSeparately.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.jkleman.mp_relax_partners_separately" );

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		using namespace protocols::docking::membrane;


		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		QuickRelaxPartnersSeparatelyOP qr( new QuickRelaxPartnersSeparately() );
		JobDistributor::get_instance()->go(qr);
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
