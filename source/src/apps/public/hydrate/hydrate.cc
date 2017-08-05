// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/hydrate/hydrate.cc
/// @brief hydrate/SPaDES protocol
/// @author Joaquin Ambia, Jason K. Lai
/// @detailed

// Protocols Headers
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/hydrate/Hydrate.hh>

// Core Headers
#include <devel/init.hh>

int
main( int argc, char * argv [] )
{
	// initialize Rosetta
	devel::init(argc, argv);

	protocols::hydrate::HydrateOP hydrate_protocol( new protocols::hydrate::Hydrate );
	protocols::jd2::JobDistributor::get_instance()->go( hydrate_protocol );

	return 0;
}
