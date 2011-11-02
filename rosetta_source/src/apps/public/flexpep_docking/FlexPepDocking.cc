// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
//
/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
#include <protocols/flexpep_docking/FlexPepDockingProtocol.hh>

#include <core/io/raw_data/ScoreFileData.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableString.hh>
// AUTO-REMOVED #include <basic/datacache/DiagnosticData.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/util.hh>

#include <protocols/moves/Mover.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
#include <algorithm>
// AUTO-REMOVED #include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


using basic::T;
using basic::Error;
using basic::Warning;

static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

//typedef utility::pointer::owning_ptr< BaseJobDistributor< BasicJobOP > > BaseJobDistributorOP
basic::Tracer TR("pilot_apps.FlexPepDock");
//static numeric::random::RandomGenerator JDRG(32342524); // magic number copied from Job Distributor

// copied from
int distribute_jobs(protocols::moves::Mover& mover, bool random_permutation)
{
  time_t overall_start_time = time(NULL);
  utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

#ifndef USEMPI
  // Reduce read contention between processes by randomizing the order in which structures are processed
  // Do not randomize, though, if job distribution is controlled by MPI
  if( random_permutation ) {
    numeric::random::random_permutation( input_jobs, numeric::random::RG );
  }
#endif

  protocols::jobdist::BasicJobOP curr_job, prev_job;
  int curr_nstruct, num_structures_processed = 0;
  core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
  core::pose::PoseOP native_pose; // starts NULL, coords *never* modified!

  protocols::jobdist::BaseJobDistributorOP jobdist;

  // pick the job distributor type based on a flag
  using basic::options::option;
  using namespace basic::options::OptionKeys;

  // What on earth is "raw" ? Surely these are called silent files ?
  bool const is_raw = option[ out::file::raw ]();
  bool const silent_output = option[ out::file::silent ].user();
  if ( is_raw || silent_output ) {
    jobdist = new protocols::jobdist::PlainRawJobDistributor(input_jobs, ".out");
  } else {
    std::string scorefile_name;
    if ( option[ out::file::scorefile ].user() ){
      scorefile_name = option[ out::file::scorefile ]();
    } else {
      scorefile_name = "score";
    }
    jobdist = new protocols::jobdist::PlainPdbJobDistributor(input_jobs, scorefile_name);
  }

  if( option[ out::nooutput ]() ){
    jobdist->disable_output();
    jobdist->enable_ignorefinished();
  }

  std::map < std::string, core::Real > score_map;

  jobdist->startup();
  while( jobdist->next_job(curr_job, curr_nstruct) ) {
    time_t pdb_start_time = time(NULL);
    TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;
    jobdist->temp_file( curr_job->output_tag(curr_nstruct) );

    // we read each PDB just once to save on disk I/O
    if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
      input_pose = new core::pose::Pose();
      if ( option[ in::file::centroid_input ].user() ) {
	core::import_pose::centroid_pose_from_pdb( *input_pose, curr_job->input_tag() );
	native_pose = new core::pose::Pose();
	core::import_pose::centroid_pose_from_pdb( *native_pose, curr_job->native_tag() );
      } else {
	core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
	native_pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *native_pose, curr_job->native_tag() );
      }
    }
    mover.set_input_pose( input_pose );
    mover.set_native_pose( native_pose );

    // Make a modifiable copy of the pose read from disk
    core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
    the_pose->data().set(
      core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
      new basic::datacache::CacheableString(curr_job->output_tag(curr_nstruct)));

    mover.apply( *the_pose );

    prev_job = curr_job; // pointer assignment, not a copy op
    num_structures_processed += 1;
    time_t pdb_end_time = time(NULL);
    TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

    score_map = protocols::jobdist::get_score_map( *the_pose );

    if ( option[ run::timer ].user() ){
      score_map["time"] = pdb_end_time - pdb_start_time;
    }

    jobdist->score_map( score_map );
    jobdist->dump_pose_and_map( curr_job->output_tag(curr_nstruct), *the_pose );
    jobdist->begin_critical_section();
    // -scorefile overrides -nooutput (otherwise, dump_pose_and_map would have taken care of this)
    jobdist->begin_critical_section();
    if ( option[ out::nooutput ]() && option[ out::file::scorefile ].user()){
      std::string scorefile_name ( option[ out::file::scorefile ] );
      std::string tag = curr_job->output_tag(curr_nstruct);
      core::io::raw_data::ScoreFileData sfd(scorefile_name);
      sfd.write_pose( *the_pose, score_map, tag );
    }
    jobdist->end_critical_section();

  } // loop over jobs and nstructs
  jobdist->shutdown();

  time_t overall_end_time = time(NULL);
  TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
  if ( num_structures_processed == 0 )
    basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	MoverOP fpDock = new flexpep_docking::FlexPepDockingProtocol(1,true, true);

	// read native pose: (TODO: look how this should be handled in Job Distributor 2)
	//	protocols::jd2::set_native_in_mover(*fpDock);
	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose = new core::pose::Pose;
		core::chemical::ResidueTypeSetCAP rsd_set;
		if ( option[ in::file::centroid_input ].user() ) {
		  core::import_pose::centroid_pose_from_pdb( *native_pose, option[ in::file::native ]() );
		} else {
		  core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
		}
		fpDock->set_native_pose( native_pose );
	}


	// run:
	protocols::jd2::JobDistributor::get_instance()->go(fpDock);
	//	protocols::jobdist::main_plain_mover( *fpDock);
	//protocols::jobdist::universal_main(*fpDock);

	exit(0);

}
