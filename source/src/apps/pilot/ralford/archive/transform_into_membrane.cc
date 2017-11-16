// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/ralford/transform_into_membrane.cc
///
/// @brief   Transform Pose Coordinates into a new membrane
/// @details
///
/// @author  Rebecca Alford (rflford12@gmail.com)
/// @note       Last Modified: 7/22/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/ddG/MembraneDDGMover.hh>
#include <protocols/membrane/ddG/Mutation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR( "apps.pilot.ralford.transform_into_membrane" );

/// @brief Main method
int
main( int argc, char * argv [] )
{

	using namespace core::scoring;
	using namespace protocols::membrane::ddG;
	using namespace protocols::jd2;

	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Create a new energy function
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012" );

		// Load in mutants for the appropriate task
		utility::vector1< MutationOP > ompA_mutations = ompA_task();
		utility::vector1< MutationOP > ompLA_mutations = ompLA_task();
		utility::vector1< MutationOP > bacteriorhodopsin_mutations = bacteriorhodpsin_task();

		// Make ddG calculation
		MembraneDDGMoverOP ddG = new MembraneDDGMover( 0.0001, sfxn, ompLA_mutations );
		JobDistributor::get_instance()->go( ddG );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

