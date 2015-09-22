// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/ncaa_utilities/MakeRotLib.cc
/// @brief MakeRotLib application
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/make_rot_lib/MakeRotLibMover.hh>

// protocols header
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>

// devel headers
#include <devel/init.hh>

// basic headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "MakeRotLib" );

int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );

		protocols::make_rot_lib::MakeRotLibMoverOP mrlm( new protocols::make_rot_lib::MakeRotLibMover() );
		protocols::jd2::JobOutputterOP nojo( new protocols::jd2::NoOutputJobOutputter() );
		protocols::jd2::JobDistributor::get_instance()->go( mrlm, nojo );

		TR << "\n+-----------------------------------------------------------------+\n"
			<<   "|                              DONE                               |\n"
			<<   "+-----------------------------------------------------------------+" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
