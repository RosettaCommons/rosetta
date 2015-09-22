// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/flo/SecondaryMatcher.cc
///
/// @brief
/// @author Florian Richter (floric@u.washington.edu)

#include <protocols/enzdes/SecondaryMatchProtocol.hh>


#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
//#include <devel/enzdes/DesignSilentStruct.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/JobDistributors.hh>

#include <utility/file/FileName.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <utility/excn/Exceptions.hh>


static THREAD_LOCAL basic::Tracer tr( "pilotapps.flo.SecondaryMatcher" );

using namespace core;

/////////////////////////////////////////

////main function
int
main( int argc, char * argv [])
{
    try {
	using namespace protocols::jobdist;

	devel::init(argc, argv);

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();

	PlainPdbJobDistributor jobdist( input_jobs );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if( option[ out::nstruct ].value() != 1 ){
		utility_exit_with_message("Value for -nstruct defined is not 1. For this protocol it doesn't make sense to process input structures more than once, so please restart with nstruct 1");
	}

	if( option[ out::nooutput ]() ){
		jobdist.disable_output();
		jobdist.enable_ignorefinished();
	}


	protocols::enzdes::SecondaryMatchProtocolOP secmatch_protocol = new protocols::enzdes::SecondaryMatchProtocol();

	basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

	BasicJobOP curr_job, prev_job;
	int curr_nstruct=0, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist.startup();

	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		tr << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
		}

		utility::file::FileName out_name( curr_job->input_tag() );
		std::string init_outtag = out_name.base();

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new basic::datacache::CacheableString( init_outtag ) );
		secmatch_protocol->apply( *the_pose );


		//NOTE: we are not writing anything in the main function, instead at the moment dumping out of
		//the secondary matcher protocol
		//jobdist.dump_pose_and_map( curr_job->output_tag(curr_nstruct), *the_pose );


		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		//std::cerr << curr_job->output_tag(curr_nstruct) << " done." << std::endl;
		tr << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs
	jobdist.shutdown();

	time_t overall_end_time = time(NULL);
	tr << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them or to use the -overwrite flag?" << std::endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}

