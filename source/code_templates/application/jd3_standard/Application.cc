// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--app_name--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#ifdef USEMPI
#include <mpi.h>
#endif

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/moves/mover_schemas.hh>
#include <protocols/jd2/parser/ScoreFunctionLoader.hh>
#include <protocols/jd2/parser/TaskOperationLoader.hh>

#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>

--mover_path--

// utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR("--app_name--");

--new_app_options_out--


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		--new_app_options_in--

		devel::init( argc, argv );

		protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();

		protocols::jd3::JobQueenOP queen( new --class--JobQueen );
		jd->go( queen );

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
