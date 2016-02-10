// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Daniel J. Mandell

// Unit Headers
#include <devel/KinematicLooprelax/KinematicLooprelax.hh>

// Rosetta Headers
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

// C++ Headers
#include <iostream>

// DJM: debug tmp
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kinematic_relax_test.main" );

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace protocols;
	using namespace protocols::jobdist;

	using basic::T;
	using basic::Error;
	using basic::Warning;

	TR << "start loops modeling" << std::endl;

	// setup the starting pose and kinematic loop relaxer once
	core::pose::Pose input_pose;
	core::import_pose::pose_from_file( input_pose, option[ OptionKeys::loops::template_pdb ]().name() , core::import_pose::PDB_file);
	devel::KinLooprelax::KinematicLooprelax kin_looprelaxer(true);

	// setup job distributor
	utility::vector1< BasicJobOP > input_jobs;
	int const nstruct_flag = option[ OptionKeys::out::nstruct ];
	int const nstruct = std::max( 1, nstruct_flag );
	std::string out_tag = option[ OptionKeys::out::output_tag ];
	BasicJobOP job = new BasicJob(out_tag, out_tag, nstruct);
	input_jobs.push_back( job );
	PlainPdbJobDistributor< BasicJobOP > jobdist( input_jobs );
	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();

	while ( jobdist.next_job(curr_job, curr_nstruct) ) { // loop over jobs
		//time_t pdb_start_time = time(NULL);
		// std::cout << "Starting " << curr_job->output_tag(curr_nstruct) << " ...\n";
		core::pose::Pose pose = input_pose;
		kin_looprelaxer.apply( pose );
		//jobdist.dump_pose_and_map( curr_job->output_tag(curr_nstruct), pose );
		// jodist.dump won't have score vs. rmsd info we want. so do this for now:
		std::string outname=option[ OptionKeys::out::path::path ]().name()+curr_job->output_tag(curr_nstruct)+".pdb";
		std::ofstream out(outname.c_str(), std::ios::out | std::ios::binary);
		core::io::pdb::dump_pdb( pose, out );
		out << "loop_rms: " << kin_looprelaxer.get_last_loop_rmsd() << std::endl;
		out << "total_energy: " << kin_looprelaxer.get_last_total_energy() << std::endl;
		out << "chainbreak: " << kin_looprelaxer.get_last_chainbreak() << std::endl;

		if (option[ OptionKeys::out::pdb_gz ]()) {
			utility::file::gzip( outname, true );
		}
		// prepare for next iteration
		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
	}
	jobdist.shutdown();
	return 0;
}

int
main( int argc, char * argv [] )
{
	try {

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );

	}	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
  return 0;
}
