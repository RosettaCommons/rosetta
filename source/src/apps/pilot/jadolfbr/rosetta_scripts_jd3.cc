// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jadolfbr/rosetta_scripts_jd3.cc
/// @brief First implementation of a JD3 app for RosettaScripts.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/moves/mover_schemas.hh>

#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>

#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR("rosetta_scripts_jd3");


using namespace protocols::rosetta_scripts;

/// Define your JobQueen
class ParsedProtocolJobQueen : public protocols::jd3::standard::StandardJobQueen
{
public:
	ParsedProtocolJobQueen()
	{

		//This is where you will list the options you will be using.  These options are command-line options
		// that can be specified at the <common> and <job> level.

		utility::options::OptionKeyList opts;

		/// Static Functions setup the OptionKeyList:
		//
		RosettaScriptsParser::list_options_read( opts );

		add_options( opts );
		parser_ = RosettaScriptsParserOP( new RosettaScriptsParser() );
	}

	~ParsedProtocolJobQueen() {}

	virtual
	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< protocols::jd3::JobResultCOP > const &
	)
	{

		TR << "Completing larval job maturation" << std::endl;

		protocols::jd3::standard::MoverAndPoseJobOP mature_job( new protocols::jd3::standard::MoverAndPoseJob );
		core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
		mature_job->pose( pose );

		bool modified_pose;

		parser_->set_recursion_limit( *job_options );
		protocols::rosetta_scripts::ParsedProtocolOP mover_protocol = parser_->generate_mover_and_apply_to_pose(*pose, *job_options, modified_pose );

		mature_job->mover( mover_protocol );
		return mature_job;
	}

	virtual bool more_jobs_remain() { return false; }

private:
	protocols::rosetta_scripts::RosettaScriptsParserOP parser_;

};

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		//Enables quick one of running as same behavior as rosetta_scripts application.

		if ( option[ parser::info ].user() ) { // If the -parser::info option is used, just print information about the requested component(s) and exit.
			protocols::rosetta_scripts::print_information( option[ parser::info ]() );
		} else if ( option[ parser::output_schema ].user() ) {
			protocols::rosetta_scripts::save_schema( option[ parser::output_schema ] );
		} else { // If an input script has been provided, then we're not printing a template script and exiting.

			//View does not make sense for JD3 - I think more would need to be done to get that compatible.  I think watching a single process should be possible,
			// but not sure what it would take to get this to work.  For now, we skip it.

			protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();

			protocols::jd3::JobQueenOP queen( new ParsedProtocolJobQueen );
			jd->go( queen );
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
