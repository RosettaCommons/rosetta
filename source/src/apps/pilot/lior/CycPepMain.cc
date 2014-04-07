#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>//option.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <protocols/loops/Loops.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/util/disulfide_util.hh>
#include <protocols/loops/LoopRelaxMover.fwd.hh>
#include <devel/init.hh>
#include <devel/cycpep/CycPepMover.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/basic.hh>
//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#define INF 99999999
// C++ headers
//#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
//
//using core::util::T;
//using core::util::Warning;
//using core::util::Error;
//namespace ss_loop_builder {
//	core::options::IntegerOptionKey start_cys("ss_loop_builder:start");
//	core::options::IntegerOptionKey   end_cys("ss_loop_builder:end");
//	core::options::IntegerOptionKey   ncycles("ss_loop_builder:ncycles");
//}
//std::string int2string(int n){
//	      std::stringstream Num;
//	      std::string str;
//	      Num << n;
//	      str = Num.str();
//	      return str;
//}
//void disulfide(core::pose::Pose& pose) {
//	using namespace core;
//	using namespace core::options;
//	using namespace core::options::OptionKeys;
//	protocols::loops::LoopRelaxMover loop_relax_mover;
//	loop_relax_mover.remodel("perturb_kic");
//	loop_relax_mover.refine("refine_kic");
//	loop_relax_mover.copy_sidechains(false);
//	Size pepsize = pose.n_residue();
//
//}
//
//int main( int argc, char** argv ){
//
//	using namespace core;
//	using namespace core::options;
//	using namespace core::options::OptionKeys;
//	option.add( ss_loop_builder::start_cys,
//								"lower seq. cys to form the loop in pose numbering").def(-1);
//		option.add( ss_loop_builder::end_cys,"higher seq. cys to form the loop in pose numbering").def(-1);
//		option.add( ss_loop_builder::ncycles,"Num. cycles of relax/loopmodel and S-S samples").def(1);
//	devel::init(argc,argv);
//
//	core::util::Tracer TR("protocols.moves.CycPep");
//
//	pose::Pose workpose;
//	pose::Pose orgpose;
//	pose::Pose native;
//	io::pdb::pose_from_pdb(orgpose, options::start_file());
//	option[ out::file::fullatom ].value(true);
//
////	std::string native_fname = option[ OptionKeys::in::file::native ];
//
////	io::pdb::pose_from_pdb( native, native_fname);
//	protocols::CycPepMover::CycPepMover cycmover;
//	cycmover.apply(orgpose);
//
//}
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jd2/JobDistributor.hh>

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


//Auto Headers
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

using basic::T;
using namespace protocols::jobdist;
using basic::Error;
using basic::Warning;
static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

//typedef utility::pointer::owning_ptr< BaseJobDistributor< BasicJobOP > > BaseJobDistributorOP
basic::Tracer TR("pilot_apps.CycPepModeler");
//static numeric::random::RandomGenerator JDRG(32342524); // magic number copied from Job Distributor

// copied from
int distribute_jobs(protocols::moves::Mover& mover, bool random_permutation)
{
  time_t overall_start_time = time(NULL);
  utility::vector1< BasicJobOP > input_jobs = load_s_and_l();

#ifndef USEMPI
  // Reduce read contention between processes by randomizing the order in which structures are processed
  // Do not randomize, though, if job distribution is controlled by MPI
  if( random_permutation ) {
    numeric::random::random_permutation( input_jobs, numeric::random::RG );
  }
#endif

  BasicJobOP curr_job, prev_job;
  int curr_nstruct, num_structures_processed = 0;
  core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
  core::pose::PoseOP native_pose; // starts NULL, coords *never* modified!
  //fdsfd
  BaseJobDistributorOP jobdist;

  // pick the job distributor type based on a flag
  using basic::options::option;
  using namespace basic::options::OptionKeys;

  // What on earth is "raw" ? Surely these are called silent files ?
  bool const is_raw = option[ out::file::raw ]();
  bool const silent_output = option[ out::file::silent ].user();
  if ( is_raw || silent_output ) {
    jobdist = new PlainRawJobDistributor(input_jobs, ".out");
  } else {
    std::string scorefile_name;
    if ( option[ out::file::scorefile ].user() ){
      scorefile_name = option[ out::file::scorefile ]();
    } else {
      scorefile_name = "score";
    }
    jobdist = new PlainPdbJobDistributor(input_jobs, scorefile_name);
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
//      if ( option[ in::file::centroid_input ].user() ) {
//	core::io::pdb::centroid_pose_from_pdb( *input_pose, curr_job->input_tag() );
//	native_pose = new core::pose::Pose();
//	core::io::pdb::centroid_pose_from_pdb( *native_pose, curr_job->native_tag() );
//      } else {
    	core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
		native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, curr_job->native_tag() );
//      }
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

    score_map = get_score_map( *the_pose );

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
    Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {

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

	MoverOP cycpepMover = new CycPepMover::CycPepMover();

	// read native pose: (TODO: look how this should be handled in Job Distributor 2)
	core::pose::PoseOP native_pose =
	  new core::pose::Pose();
	if (option[ in::file::native ].user()) {
//	  if ( option[ in::file::centroid_input ].user() ) {
//	    core::io::pdb::centroid_pose_from_pdb( *native_pose, option[ in::file::native ]() );
//	  } else {
	    core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
//	  }
	}
	cycpepMover->set_native_pose( native_pose );

	// run:
	protocols::jd2::JobDistributor::get_instance()->go(cycpepMover);
	//	protocols::jobdist::main_plain_mover( *fpDock);
	//protocols::jobdist::universal_main(*fpDock);

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}
