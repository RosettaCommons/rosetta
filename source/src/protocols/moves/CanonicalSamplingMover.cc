#include <protocols/moves/CanonicalSamplingMover.hh>
#include <protocols/moves/SidechainMover.hh>
#include <protocols/moves/SidechainMover.fwd.hh>
#include <protocols/jd2/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIBPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/HPool.hh>
#include <protocols/moves/mc_convergence_checks/HierarchicalLevel.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/PathName.hh>

#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <fstream>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/Fmath.hh>
#ifdef USEMPI
#include <mpi.h>
#endif

OPT_1GRP_KEY(Boolean, sampling, no_detailed_balance) //option moved to test_canonical_mover, need to access this value in setup
OPT_1GRP_KEY(Integer,sampling,ntrials)
OPT_1GRP_KEY(Real,sampling,mc_kt)
OPT_1GRP_KEY(Integer, sampling, interval_pose_dump)
OPT_1GRP_KEY(Integer, sampling, interval_data_dump)
OPT_1GRP_KEY(Boolean, sampling,output_only_cluster_transitions)
OPT_1GRP_KEY(Real, sampling, transition_threshold )
OPT_2GRP_KEY(File, sampling, out, new_structures )
OPT_1GRP_KEY(Integer, sampling, max_files_per_dir )

//debug ramping up temperature to equilibrate structure?
OPT_1GRP_KEY(Boolean, sampling, ramp_temperature)
OPT_1GRP_KEY(Integer, sampling, interval_increment_temp)
OPT_1GRP_KEY(Real, sampling, starting_temp)
OPT_1GRP_KEY(Boolean, sampling, add_constraints)

//dump or save part of the structure? 
OPT_1GRP_KEY(Boolean, sampling, save_loops_only )
OPT_1GRP_KEY(Boolean, sampling, dump_loops_only )

namespace protocols{
namespace moves {

static thread_local basic::Tracer tr( "protocols.moves.CanonicalSamplingMover" );

bool protocols::moves::CanonicalSamplingMover::options_registered_( false );

void CanonicalSamplingMover::register_options() {
  using namespace protocols::moves;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::pack::task;
  if ( !options_registered_ ) {
    OPT( mc::known_structures);
    OPT( out::path::path );
    NEW_OPT( sampling::ntrials, "number of Monte Carlo trials to run", 1000);
    NEW_OPT( sampling::no_detailed_balance, "preserve detailed balance", false );
    NEW_OPT( sampling::mc_kt,"value of kT for Monte Carlo",0.6);
    NEW_OPT( sampling::interval_pose_dump,"dump a pose out every x steps",1000);
    NEW_OPT( sampling::interval_data_dump,"dump data out every x steps",100);
    NEW_OPT( sampling::output_only_cluster_transitions, "output only cluster transitions", false);
    NEW_OPT( sampling::transition_threshold, "if rmsd to known_structures larger than X, add a new structure to pool", 0.5 );
    NEW_OPT( sampling::out::new_structures, "write structures above transition_threshold to this file", "discovered_decoys.out" );
    NEW_OPT( sampling::max_files_per_dir, "distribute traj and transition files into subdirectories with max N entries", 1000 );

    //
    //debug ramping up temperature to equilibrate structure?
    NEW_OPT( sampling::ramp_temperature, "ramp up the temperature and use constraints to equilibrate structure", false);
    NEW_OPT( sampling::interval_increment_temp, "increment the temperature by 0.1 every x steps", 100000 );
    NEW_OPT( sampling::starting_temp, "increment the temperature by 0.1 every x steps", 0.1 );
    NEW_OPT( sampling::add_constraints, "add constraints during equilibration?", false);

    //dump or save part of the structure?
    NEW_OPT( sampling::save_loops_only, "save only loop conformation to pool", false );
    NEW_OPT( sampling::dump_loops_only, "dump only loop conformation in silent-files" , false );

    options_registered_ = true;
  }



}

CanonicalSamplingMover::CanonicalSamplingMover():
  Mover("CanonicalSamplingMover"),
  mc_(),
  sfxn_(),
  randmove_(new RandomMover()),
  pool_rms_(),
  interval_posedump_(100),
  interval_transitiondump_(100),
  ntrials_(1000),
  detailed_balance_(true),
  MPI_synchronize_pools_(false),
  use_hierarchical_clustering_(false),
  save_loops_only_(false),
  dump_loops_only_(false),
  output_only_cluster_transition_(false)
{
  set_defaults_from_cmdline();
}

CanonicalSamplingMover::CanonicalSamplingMover(
  core::scoring::ScoreFunctionOP sfxn,
  mc_convergence_checks::Pool_RMSD_OP ptr,
  int ntrial
):
  Mover("CanonicalSamplingMover"),
  mc_(new MonteCarlo( *sfxn, basic::options::option[ basic::options::OptionKeys::sampling::mc_kt ]() ) ),
  sfxn_(sfxn),
  randmove_(new RandomMover()),
  pool_rms_(ptr),
  interval_posedump_(1000),
  interval_transitiondump_(100),
  ntrials_(ntrial),
  detailed_balance_(true),
  MPI_synchronize_pools_(false),
  use_hierarchical_clustering_(false),
  save_loops_only_(false),
  dump_loops_only_(false),
  output_only_cluster_transition_(false)
{
  set_defaults_from_cmdline();
  runtime_assert( sfxn );
}

void CanonicalSamplingMover::set_defaults_from_cmdline() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  runtime_assert( options_registered_ );

  ntrials( option[ sampling::ntrials ] );
  detailed_balance( !option[ sampling::no_detailed_balance ] );
  output_only_cluster_transitions( option[ sampling::output_only_cluster_transitions ] );
  set_interval_pose_dump( option[ sampling::interval_pose_dump ] );
  set_interval_data_dump( option[ sampling::interval_data_dump ] );
  transition_threshold_ = option[ sampling::transition_threshold ]();
  ramp_temperature_ = option[ sampling::ramp_temperature ]();
  save_loops_only_ =  option[ sampling::save_loops_only ]();
  dump_loops_only_ =  option[ sampling::dump_loops_only ]();
		    
}

void
CanonicalSamplingMover::add_mover(
  MoverOP m,
  core::Real probability
) {
  randmove_->add_mover( m, probability);
}

//copied from src/apps/pilot/dekim/bbin.cc
core::Real
CanonicalSamplingMover::periodic_range(
	       core::Real a,
	       core::Real x
	       )
{
  using namespace ObjexxFCL;
  core::Real const halfx = 0.5f * x;
  return ( ( a >= halfx || a < -halfx ) ? mod( mod( a, x ) + ( x + halfx ), x ) - halfx : a );
}



  std::string CanonicalSamplingMover::get_ABGEO_string( core::pose::Pose & pose, protocols::loops::Loops & loop ) {
 
  std::string ABGEO_assignment = "";
  for( protocols::loops::Loops::const_iterator itr = loop.begin(); itr != loop.end(); itr++ ) {
    for( core::Size ii = itr->start(); ii <= itr->stop(); ii++ ){
      core::Real phi = pose.phi(ii);
      core::Real psi = pose.psi(ii);
      core::Real omega = pose.omega(ii);
      periodic_range( phi  , 360.0 );  //does this get applied to phi??
      periodic_range( psi  , 360.0 );
      periodic_range( omega, 360.0 );
      std::string position_assignment="";
      if ( std::abs( omega ) < 90 ) {
	position_assignment= "O";
      } else if ( phi >= 0.0 ) {
	if ( -100 < psi && psi <= 100 ) {
	  position_assignment= "G"; // alpha-L
	} else {
	  position_assignment= "E"; // E
	}
      } else {
	if ( -125 < psi && psi <= 50 ) {
	  position_assignment= "A"; // helical
	} else {
	  position_assignment= "B"; // extended
	}
      }
      ABGEO_assignment = ABGEO_assignment + position_assignment;
    }
    ABGEO_assignment = ABGEO_assignment + ",";
  }
  return ABGEO_assignment;
}

void CanonicalSamplingMover::ntrials(int ntrials) {ntrials_ = ntrials;}

void CanonicalSamplingMover::set_temp(core::Real temperature) {
  mc_ = new MonteCarlo(*sfxn_, temperature);
  temperature_ = temperature;
}

void CanonicalSamplingMover::set_interval_pose_dump(int p_interval){ interval_posedump_=p_interval;}

void CanonicalSamplingMover::set_interval_data_dump(int d_interval){ interval_transitiondump_=d_interval;}

void CanonicalSamplingMover::set_scorefunction(core::scoring::ScoreFunctionOP sfxn) {sfxn_ = sfxn;}

void CanonicalSamplingMover::detailed_balance(bool truefalse) {detailed_balance_ = truefalse;}


void CanonicalSamplingMover::use_MPI_sync_pools(bool truefalse) {MPI_synchronize_pools_ = truefalse;}

void CanonicalSamplingMover::use_MPI_bcast(bool truefalse) {MPI_bcast_ = truefalse;}

void CanonicalSamplingMover::use_hierarchical_clustering(bool truefalse) {use_hierarchical_clustering_ = truefalse;}


void CanonicalSamplingMover::set_poolrmsd(mc_convergence_checks::Pool_RMSD_OP ptr){
  pool_rms_ = ptr;
}

void CanonicalSamplingMover::output_only_cluster_transitions(bool truefalse){
  output_only_cluster_transition_ = truefalse;
}

void CanonicalSamplingMover::setup_constraints( core::pose::Pose & pose ){
  pose.remove_constraints();
  core::Real const CA_cutoff(9.0);
  core::Real const cst_tol(0.5);
  for(unsigned int itr_res_i = 1; itr_res_i <= pose.total_residue(); itr_res_i++){
    for(unsigned int itr_res_j = 1; itr_res_j <= pose.total_residue(); itr_res_j++){
      Vector const CA_i( pose.residue( itr_res_i ).xyz(" CA "));
      Vector const CA_j( pose.residue( itr_res_j ).xyz(" CA "));
      core::Real const CA_dist = ( CA_i - CA_j ).length();
      if( CA_dist < CA_cutoff ){
	pose.add_constraint(
	    new core::scoring::constraints::AtomPairConstraint(
	       core::id::AtomID(pose.residue(itr_res_i).atom_index(" CA "),itr_res_i),
	       core::id::AtomID(pose.residue(itr_res_j).atom_index(" CA "),itr_res_j),
	       new core::scoring::constraints::HarmonicFunc( CA_dist, cst_tol )
	    )
	);
      }
    }
  }
}

std::string jobname_dirhash( std::string const& dir_prefix, core::Size nr_dirs ) {
  //this should split all files over <nr_dirs> directories:
  // to do this I simply remove the last 5 chars from jobname (which are probably the nstruct _0001)
  // then I add all ascii values of each character and take it modulo nr_dirs
//   core::Size sum( 0 );
//   const char* str = jobname.c_str();
//   for ( core::Size i=1; i<=jobname.size()-5; i++ ) {
//     sum+=(int) *(str++);
//   }
  core::Size job_id = protocols::jd2::JobDistributor::get_instance()->current_job_id();
  utility::file::PathName output_path = basic::options::option[ basic::options::OptionKeys::out::path::path ];
  std::string dir_name(
		       output_path.name() +"/" +
		       dir_prefix +"/" +
		       ObjexxFCL::lead_zero_string_of( job_id % nr_dirs, 4 ) +"/"
		       );
  utility::file::create_directory_recursive( dir_name );
  return dir_name;
}

void CanonicalSamplingMover::dump_decoy_or_score(
   std::ostream& os,
   core::pose::Pose const& pose,
   core::Size i_trial,
   std::string const& jobname,
   loops::Loops const& loop_to_dump,
   bool score_only /*default fasle*/
) {

   //write to silent-struct
  PROF_START( basic::CANONICALMOVER_WRITE_TO_FILE );
  core::io::silent::SilentStructOP ss;
  if ( score_only ) {
    ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "score" );
  } else {
    ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
  }

  using namespace ObjexxFCL;
  //allow two easy obtainable skip intervals: 10 and 100
  core::Size itrial_100 = i_trial / interval_posedump_ / 100;
  core::Size rest100 =  ( i_trial / interval_posedump_ ) % 100;
  core::Size rest100_10 = rest100 / 10;
  core::Size rest100_rest10 = rest100 % 10;
  core::Size score_dump100 =  ( i_trial / interval_transitiondump_ ) % interval_posedump_;
  std::string score_dump_str("");
  if ( score_only ) {
    score_dump_str="."+lead_zero_string_of( score_dump100, 3 );
  }

  ///construct NAME_00001.0.0
  ///          NAME_00001.0.1
  //               ...
  ///          NAME_00001.1.0
  ///   for easy extraction of 10th and 100th parts of structures from silent-file
  ///   10th grep '\..\.0'
  ///   100th grep '\.0\.0'

  std::string decoy_tag = jobname + "_"
    + lead_zero_string_of( itrial_100, 8 ) + "."
    + lead_zero_string_of( rest100_10, 1 ) + "."
    + lead_zero_string_of( rest100_rest10, 1)
    + score_dump_str;

  if( dump_loops_only_ /* save only loop conformations */ 
      && loop_to_dump.num_loop() > 0 ) {
    // make pose with just loop coordinates

    for( loops::Loops::const_iterator itr = loop_to_dump.begin(); itr != loop_to_dump.end(); itr++ ) {
      core::pose::Pose looponly( pose, itr->start(), itr->stop() );
      ss->fill_struct(looponly, decoy_tag);
      //looponly.copy_segment(itr->size(),pose,looponly.total_residue()+1,itr->start());
    }

  } else {
    ss->fill_struct(pose, decoy_tag);
  }
  ss->add_energy( "itrial", i_trial );
  if ( i_trial == 0 ) { //first time ?
    ss->print_header( os );
  }
  core::io::silent::SilentFileData sfd;
  tr.Debug << "about to write to silent file" << std::endl;
  sfd.write_silent_struct(*ss, os, score_only );
  PROF_STOP( basic::CANONICALMOVER_WRITE_TO_FILE );
}

void
CanonicalSamplingMover::apply(Pose & pose){
  using namespace ObjexxFCL::format;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  tr.Debug << "starting on pose with name: " << (pose.pdb_info())->name() << std::endl;

  if ( pose.total_residue() == 0 ) {
    utility_exit_with_message( "did you forget -in:file:silent ? Need to start CanonicalSamplingMover with a valid pose" );
  }

  runtime_assert( pose.total_residue() > 0 );

  PROF_START( basic::MPICANONICALSAMPLING );
  std::string jobname = protocols::jd2::current_output_name();

  runtime_assert( pool_rms_ );

  // set up loop definition if we're only sampling loop conformations. 
  // even if we're not sampling loop defs, make empty loop definition
  loops::Loops loops;
  if( option[basic::options::OptionKeys::loops::loop_file].user() ) {
    std::string loopfile = option[basic::options::OptionKeys::loops::loop_file]()[1];
    loops.read_loop_file( loopfile );
  }
  //add constraints if specified
  if( option[ constraints::cst_file ].user() ) {
    core::scoring::constraints::add_constraints_from_cmdline( pose, *sfxn_ );
  }


#ifdef USEMPI
  int n_nodes;
  MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &n_nodes ) );

  if( MPI_synchronize_pools_  && MPI_bcast_ ){
    mc_convergence_checks::MPIBPool_RMSD* mpi_pool_rms = dynamic_cast <mc_convergence_checks::MPIBPool_RMSD *> (&(*pool_rms_));
    mpi_pool_rms->set_discovered_out( (option[ sampling::out::new_structures ]()).name() );
    if( mpi_pool_rms ){
      mpi_pool_rms->set_transition_threshold( transition_threshold_ );
      //tr.Debug << " set transition threshold to " << transition_threshold_ << std::endl;
      mpi_pool_rms->set_reserve_size( n_nodes  );
      pool_rms_ = mpi_pool_rms;
    }else{
      utility_exit_with_message("cast to MPIBPool_RMSD failed! fatal error!");
    }
  }else if( MPI_synchronize_pools_  ){
    if( !use_hierarchical_clustering_ ) {
      mc_convergence_checks::MPIPool_RMSD* mpi_pool_rms = dynamic_cast <mc_convergence_checks::MPIPool_RMSD *> (&(*pool_rms_));
      mpi_pool_rms->set_discovered_out( (option[ sampling::out::new_structures ]()).name() );
      if( mpi_pool_rms ){
	mpi_pool_rms->set_transition_threshold( transition_threshold_ );
	mpi_pool_rms->set_reserve_size( n_nodes );
	mpi_pool_rms->set_reserve_size( n_nodes );
	//tr.Debug << " set transition threshold to " << transition_threshold_ << std::endl;
	pool_rms_ = mpi_pool_rms;
      }else{
	utility_exit_with_message("cast to MPIPool_RMSD failed! fatal error!");
      }
    } else {
      mc_convergence_checks::MPIHPool_RMSD* mpi_pool_rms = dynamic_cast <mc_convergence_checks::MPIHPool_RMSD *> (&(*pool_rms_));
      mpi_pool_rms->set_discovered_out( (option[ sampling::out::new_structures ]()).name() );
      if( mpi_pool_rms ){
	mpi_pool_rms->set_transition_threshold( transition_threshold_ );
	mpi_pool_rms->set_reserve_size( n_nodes );
	mpi_pool_rms->set_reserve_size( n_nodes );
	tr.Debug << "MPIHPool: set transition threshold to " << transition_threshold_ << std::endl;
	pool_rms_ = mpi_pool_rms;
	tr.Debug << "about to begin sampling with MPIHPool " << std::endl;
      }else{
	utility_exit_with_message("cast to MPIPool_RMSD failed! fatal error!");
      }
    }

  }
#endif


  //output params
  int width(10); int precision(6);
  ///make sure that this job is not repeated again in case of restart --> output to general output file
  if ( protocols::jd2::JobDistributor::get_instance()->job_outputter() ) {
    protocols::jd2::JobDistributor* jd
      = protocols::jd2::JobDistributor::get_instance();
    jd->job_outputter()->final_pose( jd->current_job(), pose );
  }

  core::Size nr_jobs( protocols::jd2::JobDistributor::get_instance()->total_nr_jobs() );
  core::Size nr_dirs( nr_jobs / option[ sampling::max_files_per_dir ]() + 1 );

  //exceptional use of ofstream to write directly to directory in MPI mode -- short-cutting the MPI-Filebuffer
  std::ofstream transition_file;
  std::string transition_filename( jobname_dirhash( "transitions",nr_dirs ) + jobname + ".transition.dat");
  transition_file.open( transition_filename.c_str() );

  ///taking care that output stream is std::ofstream so that we write directly to File instead of rerouting via MPIFileBuf.
  /// this would overload MPIFilebuffer since we open a different file for each process...
  /// keep file open... parallel open/closing is hard on the file-system
  std::string traj_filename( jobname_dirhash( "trajectories", nr_dirs ) + jobname + ".traj.out" );
  std::string traj_scorefile( jobname_dirhash( "trajectories", nr_dirs ) + jobname + ".traj.sc" );
  std::ofstream traj_file( traj_filename.c_str() );
  std::ofstream traj_sc( traj_scorefile.c_str() );
  core::io::silent::SilentFileData sfd;
  sfd.set_filename( traj_filename );

  //test-- add ramp up temperature in order to equilibrate structure?
  core::Size interval_inc_temp = option[ sampling::interval_increment_temp ];
  core::Real starting_temp = option[ sampling::starting_temp ];
  core::Real ending_temp = option[ basic::options::OptionKeys::sampling::mc_kt ];
  bool constrain_structure = option[ sampling::add_constraints ];

  if( ramp_temperature_ ){
    mc_ = new MonteCarlo( *sfxn_, starting_temp );
  }else{
    mc_ = new MonteCarlo( *sfxn_, ending_temp );
  }

  /**
  if( MPI_synchronize_pools_ && use_hierarchical_clustering ) {
    core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_silent_struct_out();
    #ifdef USEMPI
    mc_convergence_checks::MPIHPool_RMSD* hpool_ptr = dynamic_cast<mc_convergence_checks::MPIHPool_RMSD * > (&(*pool_rms_));
    if( hpool_ptr ) {
      hpool_ptr->write_headers_to_hierarchy( ss );
    }
    #endif
  }
  **/

  if ( !output_only_cluster_transition_ ) {
    transition_file << "I_TRIAL SCORE RG CLUSTER RMS_TO_CLUSTER RMS_TO_START" << std::endl;
  } else {
    transition_file << "I_TRIAL STEPS_SINCE_TRANSITON SCORE RG CLUSTER RMS_TO_CLUSTER RMS_TO_START" << std::endl;
  }

  runtime_assert( mc_ ); //are we initialized ?
  mc_->reset( pose );

  //setup Rg score calculator
  core::scoring::methods::RG_Energy_Fast rge;

  //keep a copy of current pose for rms calculations
  core::pose::Pose init( pose );
  core::Real rms_to_start;

  std::string current_cluster_center="";
  Size current_cluster_first_seen( 0 );

  if( constrain_structure ){
    setup_constraints( pose );
    sfxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
  }

    //main loop
    for ( Size i_trial=0; i_trial < ntrials_; i_trial++ ) {
      tr.Debug << "now on trial " << i_trial << " score is " << ( *sfxn_ )( pose ) << std::endl;
      core::Real proposal_density_ratio( 1 );
      randmove_->apply( pose );
      tr.Debug << "finished applying move" << std::endl;
      //(*sfxn_)(pose);
      //all just to get last_proposal_density_ratio
      if ( detailed_balance_ ) {
	proposal_density_ratio = randmove_->last_proposal_density_ratio();
      }
      tr.Debug << "got proposal density ratio: " << proposal_density_ratio << std::endl;

      mc_->boltzmann( pose, randmove_->type(), proposal_density_ratio );
	tr.Debug << "applied boltzmann criterion" << std::endl;

	//mc_->show_counters(); //debug output

      //test-- add ramp up temperature in order to equilibrate structure?
      if( ramp_temperature_ &&
	  mc_->temperature() < ending_temp &&
	  ( i_trial % interval_inc_temp ) == 0 ){

	if( constrain_structure ){
	  setup_constraints( pose );
	  sfxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	}
	mc_->set_temperature( mc_->temperature() + 0.1 );
      }
      if( ramp_temperature_ &&
	  constrain_structure &&
	  sfxn_->has_nonzero_weight( core::scoring::atom_pair_constraint ) &&
	  (i_trial % (interval_inc_temp/10) == 0)
	  ){
	sfxn_->set_weight(
			  core::scoring::atom_pair_constraint,
			  ( sfxn_->get_weight(core::scoring::atom_pair_constraint) - 0.1 )
			  ); //reduce constraints as simulation progresses
	if(sfxn_->get_weight( core::scoring::atom_pair_constraint ) < 0.1){
	  sfxn_->set_weight( core::scoring::atom_pair_constraint, 0.0); //because of numeric instability
	  pose.remove_constraints();
	}
      }
      if ( (i_trial % interval_posedump_) == 0 ) {
	tr.Debug << "about to call dump_decoy_or_score" << std::endl;
	dump_decoy_or_score( traj_file, pose, i_trial, jobname,  loops, false /*not just score*/);
      }

      if ( (i_trial % interval_transitiondump_)  == 0 ) { //output current next cluster
	std::string cluster_center; core::Real rms_to_cluster;
	core::pose::Pose looponly;
	core::Size new_level_start = 0;
	tr.Debug << "slave on " << i_trial << std::endl;
	if( use_hierarchical_clustering_ ) {
	  if( !MPI_synchronize_pools_ ) {
	    mc_convergence_checks::HierarchicalLevel* hpool_ptr = dynamic_cast<mc_convergence_checks::HierarchicalLevel * > (&(*pool_rms_));
	    utility::vector1< core::Size > address(hpool_ptr->nlevels(), 0 );
	    utility::vector1< core::Real > rms_to_cluster(hpool_ptr->nlevels(), 0.0);
	    if( save_loops_only_ && loops.num_loop() > 0 ) {
	      for( loops::Loops::const_iterator itr = loops.begin(); itr != loops.end(); itr++ ) {
		looponly = core::pose::Pose( pose, itr->start(), itr->stop() );
	      }
	      hpool_ptr->evaluate( looponly, cluster_center, rms_to_cluster, address );
	    }else {
	      hpool_ptr->evaluate( pose, cluster_center, rms_to_cluster, address );
	    }
	    mc_convergence_checks::HierarchicalLevelOP level_n = hpool_ptr;
	    bool above_threshold = false;
	    for( core::Size ii = 1; ii <= rms_to_cluster.size(); ii++ ) {
	      tr.Debug << "level " << ii <<
		"radius: " << level_n->radius() <<
		" best-rmsd: " << rms_to_cluster[ ii ] <<
		" best-cluster: " << cluster_center[ ii ] <<
		" best-address: " << address[ ii ] << std::endl;
	      if( rms_to_cluster[ ii ] > level_n->radius() ) {
		above_threshold = true;
		new_level_start = ii;
	      }
	      level_n = (level_n->next_level());
	    }

	    if( above_threshold ){
	      tr << "NEW STRUCTURE ALERT" << std::endl;
	      for( core::Size ii = 1; ii <= address.size(); ii++ ) {
		tr << address[ ii ] << " ";
	      }
	      tr << std::endl;
	      std::string newtag = "new-structure-tag";
	      if( save_loops_only_ && loops.num_loop() > 0 ){ 
		hpool_ptr->add_new( looponly, newtag, address, true, new_level_start );
	      } else {
		hpool_ptr->add_new( pose, newtag, address, true, new_level_start );
	      }
	    }
	  } else { //use hierarchy and use MPI-synching
#ifdef USEMPI
	    mc_convergence_checks::MPIHPool_RMSD* hpool_ptr = dynamic_cast<mc_convergence_checks::MPIHPool_RMSD * > (&(*pool_rms_));
	    runtime_assert( hpool_ptr != 0 );
	    tr.Debug << "about to call hierarchical mpi eval-and-add" << std::endl;
	    if( save_loops_only_ && loops.num_loop() > 0 ){
	      for( loops::Loops::const_iterator itr = loops.begin(); itr != loops.end(); itr++ ) {
		looponly = core::pose::Pose( pose, itr->start(), itr->stop() );
	      }
	      hpool_ptr->evaluate_and_add( looponly, cluster_center, rms_to_cluster);
	    } else {
	      hpool_ptr->evaluate_and_add( pose,  cluster_center, rms_to_cluster );
	    }
	    //tr << "before scoring" << std::endl;
	    //(*sfxn_)(pose);
	    //tr << "after scoring" << std::endl;
	    tr << "after MPIHPool-eval: " << cluster_center << " " << rms_to_cluster << std::endl;
#endif
	  }
	} else{ //use MPIPool
	  if( save_loops_only_ && loops.num_loop() > 0 ){
	    for( loops::Loops::const_iterator itr = loops.begin(); itr != loops.end(); itr++ ) {
	      looponly = core::pose::Pose( pose, itr->start(), itr->stop() );
	      tr << "looponly pose has " << looponly.total_residue() << " should be containing segment: " << itr->start() << " " << itr->stop() << " " << " from ref: " << pose.total_residue() <<  std::endl;	      
	    }
	    tr << "check: looponly pose has " << looponly.total_residue() << std::endl;
	    pool_rms_->evaluate_and_add( looponly, cluster_center, rms_to_cluster, transition_threshold_ );
	  } else {
	    pool_rms_->evaluate_and_add( pose, cluster_center, rms_to_cluster, transition_threshold_ );
	  }
	}
	//pool_rms_->evaluate( pose, cluster_center, rms_to_cluster);
	runtime_assert( pose.total_residue() > 0 );
	dump_decoy_or_score( traj_sc, pose, i_trial, jobname, loops, true /*not just score*/ );

	if ( !output_only_cluster_transition_ ) {
	  PROF_START( basic::MPICANONICALSAMPLING );
	  //evaluates cluster and writes transition
	  rms_to_start = core::scoring::CA_rmsd( init, pose );
	  transition_file << i_trial << " "
			  << F(width,precision,mc_->temperature()) << " "
			  << F(width,precision,(*sfxn_)( pose )) << " "
			  << F(width,precision,rge.calculate_rg_score( pose )) << " "
			  << cluster_center << " "
	                  << get_ABGEO_string( pose, loops ) << " " //add in on-the-fly ABGEO assignment to save time
			  << F(width,precision,rms_to_cluster) << " "
			  << F(width,precision,rms_to_start) << " "
			  << std::endl;
	  PROF_STOP( basic::MPICANONICALSAMPLING );
	} else {// output_only_cluster_transition_
	  //check if transition occurs
	  if ( current_cluster_center.compare("") == 0 ) {
	    //first cluster seen
	    current_cluster_center = cluster_center;
	  } else if (current_cluster_center.compare( cluster_center ) != 0 ) { //new cluster
	    PROF_START( basic::MPICANONICALSAMPLING );
	    rms_to_start = core::scoring::CA_rmsd( init, pose );
	    transition_file << i_trial << " "
			    << i_trial-current_cluster_first_seen << " "
			    << F(width,precision,mc_->temperature()) << " "
			    << F(width,precision,(*sfxn_)( pose )) << " "
			    << F(width,precision,rge.calculate_rg_score( pose )) << " "
			    << cluster_center << " "
			    << get_ABGEO_string( pose, loops ) << " " //add in on-the-fly ABGEO assignment to save time
			    << F(width,precision,rms_to_cluster)  << " "
			    << F(width,precision,rms_to_start) << " "
			    << std::endl;
	    current_cluster_center = cluster_center;
	    current_cluster_first_seen = i_trial;
	    PROF_STOP( basic::MPICANONICALSAMPLING );
	  }
	}
      } //output occured
    } //for loop over trials
    //DEBUG OUTPUT

    /**
#ifdef USEMPI
    int rank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, (int*) (&rank) );
    std::ofstream debug_cl_cnters;
    std::ostringstream q;
    q << rank;
    debug_cl_cnters.open((q.str() + ".debug_cl_centers.txt").c_str());
    for(core::Size itr = 1; itr <= pool_rms_->size(); itr++){
      std::string tag = pool_rms_->get_tag(itr);
      debug_cl_cnters << tag << std::endl;
    }
#endif
    **/
    //DEBUG OUTPUT

} //apply


std::string
CanonicalSamplingMover::get_name() const {
	return "CanonicalSamplingMover";
}


} //moves
} //protocols

