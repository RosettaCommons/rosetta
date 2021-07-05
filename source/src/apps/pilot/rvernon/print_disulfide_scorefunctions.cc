// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <devel/init.hh>



#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>

#include <utility/excn/Exceptions.hh> // AUTO IWYU For Exception







int
main( int argc, char* argv [] ) {

	try {

		// options, random initialization
		devel::init( argc, argv );

		core::scoring::disulfides::FullatomDisulfidePotentialOP potent( new core::scoring::disulfides::FullatomDisulfidePotential );

		potent->print_score_functions();

		return 0;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}
