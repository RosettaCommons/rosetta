// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/murphp/inv_kin_lig_loop_design.cc
///
/// @brief
/// @author

#include <basic/options/option.hh>
#include <core/types.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/io/pdb/pose_io.hh>
#include <basic/datacache/CacheableString.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// my devel headers
#include <utility/tag/Tag.hh>
#include <devel/inv_kin_lig_loop_design/Cloner.hh>
#include <devel/inv_kin_lig_loop_design/Protocol.hh>

#include <core/kinematics/FoldTree.hh>


// option key includes

#include <basic/options/keys/murphp.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/datacache/BasicDataCache.hh>




using namespace std;
namespace {

  void open( ifstream& fin, string const& filename ) {
    fin.close();
    fin.clear();
    fin.open(filename.c_str());
    assert( fin );
  }

  void touch( string const& filename ) {
    ofstream fout;
    fout.open( filename.c_str() );
    fout << "in_progress" << endl;
    fout.close();
  }

}



int main( int argc, char ** argv ) {
	try {
  // Parses command line options and inits RNG.
  // Doesn't seem to hurt to do it again if already done once (?)
  devel::init(argc, argv);



  // Step 0 - read in parameter file
  string const loop_design_filename = basic::options::option[ basic::options::OptionKeys::murphp::inv_kin_lig_loop_design_filename ]();//pdbfiles[1];
  cout << "loop_design_filename=" << loop_design_filename << endl;
  ifstream fin;
  open(fin,loop_design_filename);
  utility::tag::TagCOP tag = utility::tag::Tag::create(fin);
  cout << tag << endl;
  assert( tag->getName() == "loop_design" );

  // Step 0.a - read in parameters
  string const infile = tag->getTag("io")->getOption<string>("infile","in.pdb");
  int const nstruct = tag->getTag("io")->getOption<int>("nstruct",1000);
  string const output_prefix = tag->getTag("io")->getOption<string>("output_prefix","xxx");

  // Step 1 - initialize pose0 and pose1

  core::pose::PoseOP pose0( new core::pose::Pose );
  core::import_pose::pose_from_pdb( *pose0, infile );


  // Step 2 - start up job distributor

  utility::vector1< protocols::jobdist::BasicJobOP > input_jobs;
  input_jobs.push_back( protocols::jobdist::BasicJobOP( new protocols::jobdist::BasicJob( output_prefix, "yyy",  nstruct) ) );

  protocols::jobdist::PlainPdbJobDistributor jobdist( input_jobs );

  protocols::jobdist::BasicJobOP curr_job;
  int curr_nstruct;

  jobdist.startup();
  while( jobdist.next_job(curr_job, curr_nstruct) ) {

    time_t pdb_start_time = time(NULL);

    cout << "Starting " << curr_job->output_tag(curr_nstruct) << " ...\n";

    touch(curr_job->output_tag(curr_nstruct) + ".pdb");

    // NB: touching the output file...

    // Make a modifiable copy of the pose read from disk
    //core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );

    cout << "inv_kin_lig_loop_design::instantiating cloner" << endl;
    devel::inv_kin_lig_loop_design::Cloner cloner(tag,pose0);

    //size_t n = pose0->n_residue();  // unused ~Labonte
    assert( pose0->n_residue() == pose0->pdb_info()->nres() );

    cout << "inv_kin_lig_loop_design::cloning pose0" << endl;
    core::pose::PoseOP pose1 = cloner.clone();
    pose1->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new basic::datacache::CacheableString(curr_job->output_tag(curr_nstruct)));

    cout << "inv_kin_lig_loop_design::dumping clone" << endl;
    core::io::pdb::dump_pdb( *pose1, "out0.pdb" );

    cout << "inv_kin_lig_loop_design::getting fold tree" << endl;
    core::kinematics::FoldTree ft = cloner.getFoldTree();
    cout << ft << endl;
    cout << "calling pose->fold_tree" << endl;
    pose1->fold_tree( ft );

    cout << "inv_kin_lig_loop_design::setting initial configuration" << endl;
    cloner.setInitialConfig();
    pose1->energies();

    cout << "inv_kin_lig_loop_design::dumping clone with initial configuration" << endl;
    core::io::pdb::dump_pdb( *pose1, "out1.pdb" );

    cout << "inv_kin_lig_loop_design::creating protocol" << endl;
    devel::inv_kin_lig_loop_design::Protocol protocol;
    vector< devel::inv_kin_lig_loop_design::Loop > loops = cloner.getLoops();
    protocol.setParams(tag,loops);

    cout << "inv_kin_lig_loop_design::applying protocol" << endl;
    protocol.apply(pose1);

    //	  // Score new structure.  Cached energies (including *residue* energies)
    //	  // must be up-to-date in order to get sensible output.  If you remove these
    //	  // lines, you *must* insert equivalent logic at the end of all apply() methods
    //	  // (or at least for all movers that might be passed to this function).
    //	  (*scorefxn)( *the_pose );
    //	  scorefxn->accumulate_residue_total_energies( *the_pose );

    // Add to silent file
    jobdist.dump_pose_and_map( curr_job->output_tag(curr_nstruct), *pose1 );

    time_t pdb_end_time = time(NULL);
    cout << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds.\n";
  } // loop over jobs and nstructs
  jobdist.shutdown();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
  return 0;
}
