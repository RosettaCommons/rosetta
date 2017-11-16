// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/jkleman/range_relax.cc
/// @brief      Relaxes the protein by relaxing in ranges
/// @details Relaxes a protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues);
///    For the membrane version, use MPRangeRelax which runs this
///    Mover underneath
/// @author     JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/relax/RangeRelaxMover.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.jkleman.range_relax" );

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::jd2;
		using namespace protocols::relax;

		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		RangeRelaxMoverOP relax( new RangeRelaxMover() );
		JobDistributor::get_instance()->go( relax );

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
