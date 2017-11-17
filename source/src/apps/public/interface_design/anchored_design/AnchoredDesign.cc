// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/interface_design/anchored_design/AnchoredDesign.cc
/// @brief Anchored Design protocol
/// @author Steven Lewis

// Unit Headers
#include <protocols/anchored_design/AnchorMovers.hh>

// Project Headers

#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

// Numeric headers

// C++ headers

using basic::Error;
using basic::Warning;

//replaces cout
static basic::Tracer TR( "apps.public.interface_design.anchored_design.AnchoredDesign" );

/// @brief main method for anchored design.  most activity in the movers.
int
main( int argc, char* argv[] )
{
	try {
		devel::init( argc, argv );
		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new protocols::anchored_design::AnchoredDesignMover() ) );
		TR << "************************d**o**n**e**************************************" << std::endl;
	}
catch (utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
}
