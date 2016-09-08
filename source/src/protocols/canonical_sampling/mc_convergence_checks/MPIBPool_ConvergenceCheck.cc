// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifdef USEMPI
#include <protocols/canonical_sampling/mc_convergence_checks/MPIBPool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <core/pose/Pose.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <basic/Tracer.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/jd2/util.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>
#include <basic/prof.hh>

#include <utility/io/mpistream.hh>
#include <utility/exit.hh>


//SilentFileStuff
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

//option stuff
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include  <mpi.h>


namespace protocols {
namespace canonical_sampling{
namespace mc_convergence_checks {

static THREAD_LOCAL basic::Tracer tr( "MPIBPool_ConvergenceCheck" );


  core::Size MPIBPool_RMSD::master_node_;
  core::Size MPIBPool_RMSD::pool_master_node_;

using namespace ObjexxFCL;
using namespace core;
using namespace utility::io::mpi_stream;
using namespace core::io::silent;


typedef FArray2P<double> FArray2P_double;
typedef FArray2D<double> FArray2D_double;
typedef FArray3D<double> FArray3D_double;

MPI_Comm   protocols::canonical_sampling::mc_convergence_checks::MPIBPool_RMSD::MPI_COMM_POOL;

MPIBPool_RMSD::MPIBPool_RMSD( std::string const& silent_file ):
  Pool_RMSD( silent_file ),
  workers_finished_( 0 ),
  nodes_finished_( 0 ),
  pool_size_( 0 ),
  new_structures_( 0 ),
  rank_( 0 ),
  pool_rank_( 0 ),
  npes_( 0 ),
  transition_threshold_(-1),
  new_decoys_out_( "discovered_decoys.out" ),
  tracer_visible_(false),
  transfer_buf_()
{
  initialize();
  if( tracer_visible_ ){
    tr.Debug << "finished initializing!" <<std::endl;
    tr.Debug << "checking: rank " << rank_ << " has " << Pool_RMSD::size() << " structures in its pool " << std::endl;
  }
}

  /**
void
MPIBPool_RMSD::register_options(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  if( !options_registered_ ){
    //  NEW_OPT( bsampling::out::new_structures, "write structures above transition_threshold to this file", "discovered_decoys.out" );
    options_registered_ = true;
  }
}
  **/
  /**
void
MPIBPool_RMSD::set_defaults_from_cmdline(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  //runtime_assert( options_registered_ );
  //new_decoys_out_ =  option[ bsampling::out::new_structures ];
}
  **/

void
MPIBPool_RMSD::set_discovered_out( std::string const& newout ){
  new_decoys_out_ = newout;
}

std::string const&
MPIBPool_RMSD::get_discovered_out(){
  return new_decoys_out_;
}

bool MPIBPool_RMSD::is_active_node(){
  //MPI_Comm_rank( MPI_COMM_POOL, (int*) (&rank_) );
  if( tracer_visible_ ){
    tr.Debug << "checking if " << pool_rank_ << " is an active node " << !nodes_finished_[ ( pool_rank_ + 1 ) ]<< std::endl;
  }
  return !nodes_finished_[ ( pool_rank_ + 1 ) ];
}

void MPIBPool_RMSD::initialize(){

  PROF_START( basic::CHECK_COMM_SIZE );
  MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank_ ) );
  MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &npes_ ) );
  PROF_STOP( basic::CHECK_COMM_SIZE );
  PROF_START( basic::INITIALIZE );
  pool_rank_ = rank_;
  pool_npes_ = npes_;

  //assume master-node is always first in active_nodes list
  pool_master_node_ = 0;

  int new_size = npes_ - master_node_;
  transfer_buf_.nresidues_ =  Pool_RMSD::natom();
  transfer_buf_.set_size( new_size );
  //set_defaults_from_cmdline();
  if ( rank_ == master_node_ ) {
    pool_size_ = Pool_RMSD::size();
    Pool_RMSD::clear();
    assert( Pool_RMSD::size() == 0 );
  } else {
    pool_size_ = Pool_RMSD::size();
  }
  PROF_STOP( basic::INITIALIZE );
  //create new MPI_COMM_WORLD based on sub-set of nodes
  PROF_START( basic::MPICOMMCREATION );
  int index = 0;
  for(int ii = master_node_; ii < static_cast<int>(npes_); ii++){
    (transfer_buf_.int_buf1_)[ index++ ] = ii;
  }

  //initialize all num_slave dependent buffers for MPI transfers
  MPI_Group pool_group, all;
  int returnval;
  //int world_rank;
  //int new_rank;

  returnval = MPI_Comm_group( MPI_COMM_WORLD, &all);
  if ( returnval != MPI_SUCCESS ) {
    utility_exit_with_message("failed in creating a new communicator!");
  }

  returnval = MPI_Group_incl( all, (new_size), transfer_buf_.int_buf1_, &pool_group );
  if ( returnval != MPI_SUCCESS ) {
    utility_exit_with_message("failed in creating a new communicator!");
  }

  returnval = MPI_Comm_create( MPI_COMM_WORLD, pool_group, &MPI_COMM_POOL );
  if ( returnval != MPI_SUCCESS ) {
    utility_exit_with_message("failed in creating a new communicator!");
  }

  update_ranks( transfer_buf_.int_buf1_, (new_size) );
  PROF_STOP( basic::MPICOMMCREATION );

  tracer_visible_ = tr.visible();
}


void MPIBPool_RMSD::reformat( core::pose::Pose const& pose, std::string & new_tag ){
  PROF_START( basic::FARRAY_MANIPULATION );
  protocols::toolbox::fill_CA_coords( pose, transfer_buf_.temp_coord_for_evaluation_ ); //index starts at 1.
  PROF_STOP( basic::FARRAY_MANIPULATION );
  //assign new tag based on olli's scheme
  assign_tag( new_tag, 0 );
}

  void MPIBPool_RMSD::assign_tag( std::string& new_tag, core::Size optional_id_num ){
  //std::string jobname = protocols::jd2::current_output_name();
  if( rank_ == master_node_ ){
    if( tracer_visible_ ){
      tr.Debug << "assigning a tag with value " << lead_zero_string_of( ( pool_size_ + new_structures_ ) , 8 ) << std::endl;
    }
    new_tag = "new."+lead_zero_string_of( pool_size_ + new_structures_, 8 ); //+".0"+"_"+jobname
  }else{
    if( optional_id_num == 0 ){
      if( tracer_visible_ ){
	tr.Debug << "assigning a tag with value " << lead_zero_string_of( Pool_RMSD::size(), 8 ) << std::endl;
      }
      new_tag = "new."+lead_zero_string_of( Pool_RMSD::size(), 8 ); //+".0"+"_"+jobname
    }else{
      if( tracer_visible_ ){
	tr.Debug << "assigning a tag with value " << lead_zero_string_of( optional_id_num, 8 ) << std::endl;
      }
      new_tag = "new."+lead_zero_string_of( optional_id_num, 8 ); //+".0"+"_"+jobname

    }
  }
}

void MPIBPool_RMSD::increment_pool_size( core::Size num_to_add ){
  new_structures_ += num_to_add;
}


void MPIBPool_RMSD::broadcast_newest_coords( int num_to_send ){
  if ( num_to_send == 0 ) return;

  PROF_START( basic::MPIBARRIER );
  MPI_Barrier( MPI_COMM_POOL );
  PROF_STOP( basic::MPIBARRIER );

  if ( rank_ == master_node_ ) {

    assert( (int)(new_structures_) >= num_to_send );
    //core::Size current_size = new_structures_;
    if( tracer_visible_ ) {
      tr.Debug << "broadcasting " << num_to_send << " structures " << std::endl;
      for( core::Size ii = 0; static_cast<int>(ii) < num_to_send; ii++) {
	tr.Debug << " sending coordinates starting at index " << transfer_buf_.int_buf1_[ ii ] << std::endl;
      }
    }
    //    PROF_START( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    //    MPI_Bcast( &num_to_send, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    //    PROF_STOP( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    //take newest coords and put back in farray_coord_ptr_

    PROF_START( basic::COPY_COORDS );
    int shifted_index = 0;
    //int_buf1_ contains starting indices of structures you wish to send
    for( core::Size ii = 0; static_cast<int>(ii) < num_to_send; ii++ ) {
      // delete coordinates that are not being saved by shifting coordinates over to the left
      for( core::Size jj = 0; jj < (3 * transfer_buf_.nresidues_ ); jj++ ) {
	transfer_buf_.farray_coord_ptr_[ shifted_index++ ] = transfer_buf_.farray_coord_ptr_[ jj + transfer_buf_.int_buf1_[ ii ] ];
      }
    }
    PROF_STOP( basic::COPY_COORDS );
  } //master copy coordinates

  PROF_START( basic::MPI_MASTER_BCAST_COORDS );
  MPI_Bcast(
	    transfer_buf_.farray_coord_ptr_,
	    ( num_to_send * transfer_buf_.nresidues_ * 3 ),
	    MPI_DOUBLE,
	    pool_master_node_,
	    MPI_COMM_POOL
  );
  PROF_STOP( basic::MPI_MASTER_BCAST_COORDS );

  if ( rank_ != master_node_ ) { //slave

    int num_to_receive = num_to_send;
    //    PROF_START( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    //    MPI_Bcast( &num_to_receive, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    //    PROF_STOP( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    if ( tracer_visible_ ) {
      tr.Debug << "receiving " << num_to_receive << " structures " << std::endl;
    }
    //if( tracer_visible_ ){
    //tr.Debug << "outputting the received coordinates " << (num_to_receive*transfer_buf_.nresidues_*3) << std::endl;
    //for(core::Size ii = 0; ii < (num_to_receive*transfer_buf_.nresidues_*3); ii++) {
    //tr.Debug << transfer_buf_.farray_coord_ptr_[ ii ] << " ";
    //}
    //tr.Debug << std::endl;
    //}
    PROF_START( basic::COPY_COORDS );
    increment_pool_size( num_to_receive );
    core::Size index = 0;
    while( num_to_receive > 0 ){
      std::string tag_to_get;
      assign_tag( tag_to_get, 0 );
      if( tracer_visible_ ){
	tr.Debug << "assigned a new tag to new structure " << tag_to_get << " " << num_to_receive << std::endl;
      }
      array_to_farray( index );
      if( tracer_visible_ ){
	tr.Debug << "successfully assign array to farray " << std::endl;
      }
      Pool_RMSD::add( (transfer_buf_.temp_coord_for_evaluation_), transfer_buf_.nresidues_, tag_to_get );
      if( tracer_visible_ ){
	tr.Debug << "successfully added structure to pool " << std::endl;
      }
      num_to_receive--;
      index += (transfer_buf_.nresidues_ * 3);
      if( tracer_visible_ ){
	tr.Debug << "next will be accessing index " << index << std::endl;
      }
    }
    PROF_STOP( basic::COPY_COORDS );
    tr.Debug << "tabulated pool size is " << new_structures_ << " and real size is " << Pool_RMSD::size() << std::endl;
  } //slave post procession
}


void
MPIBPool_RMSD::master_gather_new_coords(){
  runtime_assert( rank_ == master_node_ );

  if( tracer_visible_ ){
    tr.Debug << "expecting " << pool_npes_ << " updates" << std::endl;
  }

  PROF_START( basic::MPI_GATHER_BARRIER );
  MPI_Barrier( MPI_COMM_POOL );
  PROF_STOP( basic::MPI_GATHER_BARRIER );

  //find out wether node has structure to report or finished trajectory
  //return in size_per_coord: -1 finished, 0 no structure, nresidues_ a structure to report
  int structures_to_report = 0;
  PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
  MPI_Gather( &structures_to_report, 1, MPI_INT, transfer_buf_.size_per_coords_, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
  PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );

  core::Size max_coord_size = 0;
  //core::Size num_added = 0;
  transfer_buf_.size_ = 0;

  for ( core::Size ii = 0; ii < ( pool_npes_ ); ii++){
    if ( tracer_visible_ ) {
      tr.Debug << "rank " << ii << " reports a size of " << transfer_buf_.size_per_coords_[ ii ] << std::endl;
      tr.Debug << "max_coord_size now has a value of " << max_coord_size << std::endl;
    }
    //read out gathered information - finished?, structure? etc.
    if ( transfer_buf_.size_per_coords_[ ii ] < 0 ){ //finished ?

      workers_finished_++;
      nodes_finished_[ ( ii + 1) ] = true;

      if ( tracer_visible_ ){
	tr.Debug << "tabulating another trajectory finished! " << std::endl;
      }
      transfer_buf_.size_per_coords_[ ii ] = 0;
      transfer_buf_.memory_offset_[ ii ] = 0;

    } else if ( transfer_buf_.size_per_coords_[ ii ] > 0 ) { //not finished --- structure to report ?
      transfer_buf_.winning_ranks_[ transfer_buf_.size_ ] = ii; //save winning rank
      transfer_buf_.memory_offset_[ ii ] = transfer_buf_.size_ * transfer_buf_.size_per_coords_[ ii ];
      transfer_buf_.size_++;
    }
    max_coord_size += transfer_buf_.size_per_coords_[ ii ];
  }

  if( tracer_visible_ ){
    tr.Debug << "checking the contents of nodes_finished: ";
    for(core::Size ii = 0; ii < pool_npes_; ii++){
      tr.Debug << nodes_finished_[ (ii+1) ] << " ";
    }
    tr.Debug << std::endl;
    tr.Debug << "about to receive the new coordinates " << max_coord_size << std::endl;
  }


  ///now receive all coordinates of announced structures
  double tmp;
  PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  MPI_Gatherv( &tmp, 0, MPI_DOUBLE, transfer_buf_.farray_coord_ptr_, transfer_buf_.size_per_coords_, transfer_buf_.memory_offset_, MPI_DOUBLE, (pool_master_node_), MPI_COMM_POOL );
  PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  if( tracer_visible_ ){
    tr.Debug << "expecting a total size of max: " << max_coord_size  << std::endl;
    tr.Debug << "received all-coordinates! ";
    for(core::Size ii = 0; ii < max_coord_size; ii++){
      tr.Debug << transfer_buf_.farray_coord_ptr_[ ii ] << " ";
    }
    tr.Debug << std::endl;
    tr.Debug << "finished receiving all coordinated" << std::endl;
  }
}


void
MPIBPool_RMSD::slave_gather_new_coords(){

  int size_to_report = transfer_buf_.temp_coord_for_evaluation_.u1() * transfer_buf_.temp_coord_for_evaluation_.u2();

  //double* array_xyz = new double[ size_to_report ];
  farray_to_array( 0 ); //prof statement in function call

  int dummy;

  if(tracer_visible_){
    tr.Debug << " calling gather to report slave size of coords" << std::endl;
  }

  PROF_START( basic::MPI_GATHER_BARRIER );
  MPI_Barrier( MPI_COMM_POOL );
  PROF_STOP( basic::MPI_GATHER_BARRIER );


  PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
  MPI_Gather( &size_to_report, 1, MPI_INT, transfer_buf_.size_per_coords_, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
  PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );

  if( tracer_visible_ ){
    tr.Debug << " slave: finished calling gather! sending coordinate of size " << size_to_report << std::endl;
  }
  PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  MPI_Gatherv(
	      transfer_buf_.farray_coord_ptr_,
	      size_to_report,
	      MPI_DOUBLE,
	      &dummy,
	      &dummy,
	      &dummy,
	      MPI_DOUBLE,
	      (pool_master_node_),
	      MPI_COMM_POOL
	      );
  PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );
}

void
MPIBPool_RMSD::slave_report_no_new_coords(){
  int num_to_report = 0;
  int dummy;
  double empty_coords;


  PROF_START( basic::MPI_GATHER_BARRIER );
  MPI_Barrier( MPI_COMM_POOL );
  PROF_STOP( basic::MPI_GATHER_BARRIER );

  if( tracer_visible_ ){
    tr.Debug << " slave: reporting no new coordinates" << std::endl;
  }
  PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
  MPI_Gather( &num_to_report, 1, MPI_INT, 0, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL);
  PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );
  PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  MPI_Gatherv( &empty_coords, 0, MPI_DOUBLE, &dummy, &dummy, &dummy, MPI_DOUBLE, (pool_master_node_), MPI_COMM_POOL);
  PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );
}

void MPIBPool_RMSD::farray_to_array( core::Size index, core::Size ){
  PROF_START( basic::FARRAY_MANIPULATION );
  for (int ii = 1; ii <= transfer_buf_.temp_coord_for_evaluation_.u1(); ii++ ) {
    for (int jj = 1; jj <= transfer_buf_.temp_coord_for_evaluation_.u2(); jj++ ) {
      transfer_buf_.farray_coord_ptr_[index++] = transfer_buf_.temp_coord_for_evaluation_( ii, jj );
    }
  }
  PROF_STOP( basic::FARRAY_MANIPULATION );
}

void MPIBPool_RMSD::farray_to_array( core::Size index ){
  farray_to_array( index, 1 );
}

void MPIBPool_RMSD::array_to_farray( core::Size index, core::Size ){
  PROF_START( basic::FARRAY_MANIPULATION );
  for ( int ii = 1; ii <= transfer_buf_.temp_coord_for_evaluation_.u1(); ii++ ) {
    for( int jj = 1; jj <= transfer_buf_.temp_coord_for_evaluation_.u2(); jj++ ) {
      transfer_buf_.temp_coord_for_evaluation_( ii, jj ) =  transfer_buf_.farray_coord_ptr_[index++];
      //tr.Debug << "outputting useful debug info: " << transfer_buf_.farray_coord_ptr_[ index - 1 ] << " " << transfer_buf_.temp_coord_for_evaluation_( ii, jj ) << std::endl;
    }
  }
  PROF_STOP( basic::FARRAY_MANIPULATION );
}

void MPIBPool_RMSD::array_to_farray( core::Size index ){
  array_to_farray( index, 1 );
}

  /**
void MPIBPool_RMSD::farray_to_array( FArray2D<double> const& farray_xyz, double xyz[] ){
  int index = 0;
  if( tracer_visible_ ){
    tr.Debug << "converting farray to array: u1: " << farray_xyz.u1() << " u2: " << farray_xyz.u2() << " " << std::endl;
  }
  for( int i = 1; i <= farray_xyz.u1(); i++ ){
    for( int j = 1; j <= farray_xyz.u2(); j++ ){
      xyz[ index++ ] = farray_xyz( i, j );
      //tr.Debug << farray_xyz( i, j ) << " ";
    }
  }
  //tr.Debug << std::endl;
}
  **/ //not needed anymore. just access memory directly

  /**
void MPIBPool_RMSD::array_to_farray( FArray2D<double>& farray_xyz, double xyz[] ){

  assert( transfer_buf_.nresidues_ > 0 );
  //farray_xyz.redimension( 3, nresidues_, 0.0 );
  tr.Debug << "converting array to farray dimensions: " << farray_xyz.u1() << " " << farray_xyz.u2() << " " << std::endl;
  int index = 0;
  for(core::Size i = 1; i <= 3; i++ ){
    for(core::Size j = 1; j <= transfer_buf_.nresidues_; j++ ){
      farray_xyz( i, j ) = xyz[ index++ ];
      //tr.Debug << farray_xyz( i, j ) << " ";
    }
  }
  //  tr.Debug << std::endl;
}
  **/


void MPIBPool_RMSD::set_transition_threshold( core::Real threshold ){
  transition_threshold_ = threshold;
}

bool MPIBPool_RMSD::workers_finished(){
  if ( workers_finished_ < ( npes_ - master_node_ - 1 ) ){
    if  ( tracer_visible_ ){
      tr.Debug << "num trajectories finished: " <<
	workers_finished_ << " needed: " <<
	( npes_ - master_node_ - 1 ) << std::endl;
    }
    return false;
  }else{
    if( tracer_visible_ ){
      tr.Debug << "FINISHED! num trajectories finished: " <<
	workers_finished_ << " needed: " <<
	( npes_ - master_node_ - 1 ) << std::endl;
    }
    return true;
  }

}

void MPIBPool_RMSD::finalize(){
  if( rank_ != master_node_ ){
    if( tracer_visible_ ){
      tr.Debug << "sending broadcast finalized message to master " << std::endl;
    }
    int size_to_report = -1;
    int empty_size = 0;
    double empty_coords;

    PROF_START( basic::MPI_GATHER_BARRIER );
    MPI_Barrier( MPI_COMM_POOL );
    PROF_STOP( basic::MPI_GATHER_BARRIER );

    PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
    MPI_Gather( &size_to_report, 1, MPI_INT, &size_to_report, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );

    PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
    MPI_Gatherv( &empty_coords, 0, MPI_DOUBLE, &empty_size, &empty_size, &empty_size, MPI_DOUBLE, (pool_master_node_), MPI_COMM_POOL );
    PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );

    int num_poses_added = 0;
    int new_size = 0;

    PROF_START( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    MPI_Bcast( &num_poses_added, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    PROF_STOP( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );

    PROF_START( basic::MPI_MASTER_BCAST_WINNING_RANKS );
    MPI_Bcast( transfer_buf_.int_buf1_, num_poses_added, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    PROF_STOP( basic::MPI_MASTER_BCAST_WINNING_RANKS );

    broadcast_newest_coords( num_poses_added );

    PROF_START( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );
    MPI_Bcast( &new_size, 1, MPI_INT, pool_master_node_, MPI_COMM_POOL );
    PROF_STOP( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );

    assert( new_size < static_cast<int>(pool_npes_) );
    if( tracer_visible_ ){
      tr.Debug << "new size is " << new_size << " current size: " << pool_npes_ << std::endl;
    }
    PROF_START( basic::COMM_REDUCE_SIZE );
    MPI_Bcast( transfer_buf_.int_buf1_, new_size, MPI_INT, pool_master_node_, MPI_COMM_POOL );
    PROF_STOP( basic::COMM_REDUCE_SIZE );

    if( tracer_visible_ ){
      tr.Debug << "creating new communicator from ranks: ";
      for(int ii = 0; ii < new_size; ii++){
	tr.Debug << transfer_buf_.int_buf1_[ ii ] << " ";
      }
      tr.Debug << std::endl;
    }

    create_comm( transfer_buf_.int_buf1_, new_size );

    //DEBUG
    /**
    std::string debug_posedump = ".debug_posedump.out";
    std::ofstream debug;
    std::ostringstream q;
    q << rank_;
    debug.open((q.str() + debug_posedump).c_str(),std::fstream::app);

    for( core::Size ii = 1; ii <= Pool_RMSD::size(); ii++ ) {
      FArray2D<double> tmp;
      std::string tag = Pool_RMSD::get_tag( ii );
      Pool_RMSD::get( ii, tmp );
      //write to file
      debug << tag << "  BEGIN " << std::endl;
      for( core::Size f_index_i = 1; f_index_i <= tmp.u2(); f_index_i++ ){
	for( core::Size f_index_j = 1; f_index_j <= tmp.u1(); f_index_j++ ){
	  debug << tmp( f_index_j, f_index_i ) << " ";
	}
	debug << std::endl;
      }
      debug << tag << " END " << std::endl;
    }
    **/
    //DEBUG

  }
}


  void MPIBPool_RMSD::create_comm( int ranks_to_include[], int new_size ){
    int returnval;
    MPI_Group new_pool_group, old_pool_group;
    MPI_Comm dup_pool_comm;
    //tr.Debug << "creating a duplicate communicator from ranks: " << std::endl;
    bool is_active_node = false;
    nodes_finished_.resize( new_size );
    for(int ii = 0; ii < new_size; ii++ ){
      nodes_finished_[ ii + 1 ] = false;
      if( (int)(pool_rank_) == ranks_to_include[ ii ]){
	tr.Debug << "this rank " << pool_rank_ << " is designated an active node" << std::endl;
	is_active_node = true;
      }
      //tr.Debug << ranks_to_include[ ii ] << " ";
    }
    //tr.Debug << std::endl;
    PROF_START( basic::MPICOMMCREATION );
    MPI_Comm_dup( MPI_COMM_POOL, &dup_pool_comm );
    returnval = MPI_Comm_group( dup_pool_comm, &old_pool_group );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved in MPIBPool_RSD::create_comm()!" );
    //tr.Debug << "created comm-group based on old pool" << std::endl;
    returnval = MPI_Group_incl( old_pool_group, (new_size), ranks_to_include, &new_pool_group );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved in MPIBPool_RSD::create_comm()!" );
    //tr.Debug << " created new group based on trajs that are still active " << std::endl;
    returnval = MPI_Comm_create( dup_pool_comm, new_pool_group, &MPI_COMM_POOL );
    runtime_assert_string_msg(returnval == MPI_SUCCESS, "MPI_SUCCESS not achieved in MPIBPool_RSD::create_comm()!"  );
    //tr.Debug << "created new comm based on this new group " << std::endl;
    if( is_active_node ){
      update_ranks( ranks_to_include, new_size );
      MPI_Comm_size(MPI_COMM_POOL, ( int* )( &new_size ) );
      ///tr.Debug << "successfully created new com. checking size: " << new_size << std::endl;
      //transfer_buf_.set_size( new_size );
      //runtime_assert( transfer_buf_.size_ > pool_npes_ );
      if( tracer_visible_ ){
	tr.Debug << "new size of pool is " << new_size << std::endl;
      }
    }
    PROF_STOP( basic::MPICOMMCREATION );
  }


/// @detail update the rank of the pool to the MPI_COMM_POOL relative rank
void MPIBPool_RMSD::update_ranks( int const active_nodes[], int new_size ){
  //bool is_active = false;

  //debug output
  if( tracer_visible_ ){
    tr.Debug << "listing active ranks: ";
    for( int ii = 0; ii < new_size; ii++ ){
      tr.Debug << active_nodes[ ii ] << " ";
    }
    tr.Debug << std::endl;
  }

  transfer_buf_.set_size( new_size );
  //figure out if this rank take parts in the POOL
  nodes_finished_.resize( new_size, false );
  for ( int ii = 0; ii < new_size; ii++ ){
    //nodes_finished_[ ii + 1 ] = false;
    if ( active_nodes[ ii ] == (int)(pool_rank_) ) {
      if( tracer_visible_ ){
	tr.Debug << "this node " << pool_rank_ << " is still active " << std::endl;
      }
      PROF_START( basic::CHECK_COMM_SIZE );
      MPI_Comm_size(MPI_COMM_POOL, ( int* )(&pool_npes_) );
      MPI_Comm_rank(MPI_COMM_POOL, ( int* )(&pool_rank_) );
      PROF_STOP( basic::CHECK_COMM_SIZE );
      if ( tracer_visible_ ){
	tr.Debug << "master node is rank " << master_node_ << std::endl;
	tr.Debug << "new size of comm is " << npes_ << " pool-size " << pool_npes_ << " this node now has rank " << rank_ << " and pool_rank " << pool_rank_ << std::endl;
	tr.Debug << "double checking node_finished contents: \n";
	for(core::Size ii = 1; ii <= nodes_finished_.size(); ii++){
	  tr.Debug << nodes_finished_[ii] << " ";
	}
	tr.Debug << std::endl;
      }
      runtime_assert( transfer_buf_.size_ >= pool_npes_ );
      return;
    }
  }
}

TransferBuffer::TransferBuffer():
  memory_offset_(0),
  size_per_coords_(0),
  int_buf1_(0),
  winning_ranks_(0),
  farray_coord_ptr_(0),
  temp_coord_for_evaluation_(),
  coords_(),
  size_(0),
  nresidues_(0)
{}

TransferBuffer::TransferBuffer( core::Size num_slave_nodes ):
  size_( num_slave_nodes )
{
  set_size( num_slave_nodes );
}

TransferBuffer::~TransferBuffer() {
  delete [] memory_offset_;
  delete [] size_per_coords_;
  delete [] int_buf1_;
  delete [] winning_ranks_;
  delete [] farray_coord_ptr_;
}

void
TransferBuffer::set_size( int num_slave_nodes ){
  if( tr.visible() ){
    tr.Debug << "setting the size of the transfer_buf_ to " << num_slave_nodes << std::endl;
  }
  memory_offset_ = new int[ num_slave_nodes ];
  size_per_coords_ = new int[ num_slave_nodes ];
  int_buf1_ = new int[ num_slave_nodes ];
  winning_ranks_ = new int[ num_slave_nodes ];
  //runtime_assert( nresidues_ > 0 && num_slave_nodes > 0);
  coords_ = ObjexxFCL::FArray3D<double>( 3, nresidues_, num_slave_nodes, 0.0 );
  temp_coord_for_evaluation_ = ObjexxFCL::FArray2D<double>( 3, nresidues_, 0.0 );
  //farray_coord_ptr_ = coords_.get_pointer_to_data();
  farray_coord_ptr_ = new double[  3 * nresidues_ * num_slave_nodes ];
  size_ = num_slave_nodes;
}


void
MPIBPool_RMSD::add_pose_to_pool(){
  std::string tag;
  assign_tag( tag, 0 );
  if( tracer_visible_ ){
    tr.Debug << "now adding a pose with the assigned-tag: " << tag << std::endl;
  }
  Pool_RMSD::add(
		 transfer_buf_.temp_coord_for_evaluation_,
		 transfer_buf_.nresidues_,
		 tag
		 );
  increment_pool_size( 1 );
}


void MPIBPool_RMSD::master_go(){

  while( !workers_finished() ){

    //using broadcasting
    if ( tracer_visible_ ){
      tr.Debug << "about to gather coords from slaves" << std::endl;
    }

    //utility::vector1<FArray2D_double> new_poses; //becomes transfer_buf_.coords
    //utility::vector1<int> rank_of_pose_added; //becomes transfer_buf_.winning_ranks

    //figure out how many slave report structures and which workers are finished
    master_gather_new_coords();

    if( tracer_visible_ ){
      tr.Debug << "finished gathering coordinates from slave nodes checking size of transfer_buf " << transfer_buf_.size_ << std::endl;
    }

    PROF_START( basic::MPI_POOL_MASTER_THINKS );
    //figure out which structures are going into the masters pool, and which are discarded
    core::Size num_poses_added = 0;
    for ( core::Size index_new_pose = 0; index_new_pose < transfer_buf_.size_; index_new_pose++) {
      if ( num_poses_added == 0 ){ //first will always go in!
	//std::string new_tag;
	//put coordinates in temp Farray
	array_to_farray( 0 );
	add_pose_to_pool();
	transfer_buf_.int_buf1_[ num_poses_added ] = 0; //offset for first structure
	num_poses_added++;
	if( tracer_visible_ ){
	  tr.Debug << "continuing . . . " << std::endl;
	}
	continue;

      }
      //for all other structure evaluate if there is already a similar one in pool.
      if( tracer_visible_ ){
	tr.Debug << "index is now " << index_new_pose << " size is " << transfer_buf_.size_ << " performing evaluation now. num structures in pool "
		 << new_structures_ << " num added in-between "
		 << num_poses_added << " so index is " << (new_structures_ - num_poses_added + 1 ) << std::endl;
      }

      core::Real best_rmsd;
      std::string best_decoy;
      if( tracer_visible_ ){
	tr.Debug << "about to convert array to farray, index:  " << (index_new_pose*transfer_buf_.nresidues_ * 3 ) << std::endl;
      }
      array_to_farray( index_new_pose * transfer_buf_.nresidues_ * 3 );
      if( tracer_visible_ ){
	tr.Debug << "finished converting array to farray, index: " << (index_new_pose) << std::endl;
      }
      Pool_RMSD::evaluate( transfer_buf_.temp_coord_for_evaluation_,
			   best_decoy,
			   best_rmsd,
			   (new_structures_ + 1 - num_poses_added )
			   );

	    //best_rmsd = transition_threshold_ + 1;

      if( tracer_visible_ ){
	tr.Debug << "finished evaluating decoy against pool, index "
		 << (new_structures_ + 1 - num_poses_added) << " best_rms: " << best_rmsd << std::endl;
      }
      if ( best_rmsd > transition_threshold_ ) {
	if ( tracer_visible_ ) {
	  tr.Debug << "finished eval, best rms is " << best_rmsd
		   << " which is greater than threshold " << transition_threshold_ << " so adding to pool " << std::endl;
	}
	add_pose_to_pool();
	transfer_buf_.int_buf1_[ num_poses_added ] = ( index_new_pose * 3 * transfer_buf_.nresidues_ ); //starting index of saved structure in array
	if( tracer_visible_ ){
	  tr.Debug << "finished adding pose to pool, now have " << new_structures_ << std::endl;
	}
	num_poses_added++;
      }
    }

    PROF_STOP( basic::MPI_POOL_MASTER_THINKS );

    //broadcast the ranks whose structure got accepted
    PROF_START( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
    MPI_Bcast( &num_poses_added, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL ); //int* ranks_with_accepted_poses
    PROF_STOP( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );

    if( tracer_visible_ ){
      tr.Debug << "broadcasting winning ranks " << std::endl;
    }
    PROF_START( basic::MPI_MASTER_BCAST_WINNING_RANKS );
    MPI_Bcast( transfer_buf_.winning_ranks_, num_poses_added, MPI_INT, (pool_master_node_), MPI_COMM_POOL ); //int* ranks_with_accepted_poses
    if( tracer_visible_ ){
      tr.Debug << "checking the contents of int_buf, which should contain starting indices\n";
      for(core::Size ii = 0; ii < num_poses_added; ii++){
	tr.Debug << transfer_buf_.int_buf1_[ ii ] << " ";
      }
      tr.Debug << std::endl;
    }
    PROF_STOP( basic::MPI_MASTER_BCAST_WINNING_RANKS );
    broadcast_newest_coords( num_poses_added );

    //after broadcast newest coords, int_buf1_ not needed anymore and can be overwritten in the next section

    //if some trajectories finished but not others, we need to re-create a comm with appropriate size
    PROF_START( basic::CHECK_COMM_SIZE );
    unsigned int new_size = 0;
    for ( unsigned int ii = 0; ii < pool_npes_; ii++ ){
      if ( !nodes_finished_[ (ii + 1) ] ) {
	transfer_buf_.int_buf1_[ new_size ] = ii;
	new_size++;
	if( tracer_visible_ ){
	  tr.Debug << "still active node: " << transfer_buf_.int_buf1_[ new_size - 1 ] << std::endl;
	}
      }
    }
    PROF_STOP( basic::CHECK_COMM_SIZE );
    PROF_START( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );
    MPI_Bcast( &new_size, 1, MPI_INT, pool_master_node_, MPI_COMM_POOL );
    PROF_STOP( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );

    if( new_size != pool_npes_ ){
      PROF_START( basic::MPI_MASTER_BCAST_NEW_POOL_RANKS );
      MPI_Bcast( transfer_buf_.int_buf1_, new_size, MPI_INT, (pool_master_node_), MPI_COMM_POOL ); //here's where we re-create the mpi_pool_comm
      PROF_STOP( basic::MPI_MASTER_BCAST_NEW_POOL_RANKS );
      if( tracer_visible_ ){
	tr.Debug << "new size is " << new_size << " current size is " << pool_npes_ <<std::endl;
      }
      PROF_START( basic::MPICOMMCREATION );
      create_comm( transfer_buf_.int_buf1_, new_size );
      PROF_STOP( basic::MPICOMMCREATION );
    }

  }
    //check num_to_update. if non-zero, broadcast numm
  if( tracer_visible_ ){
    tr.Debug << "master node finished " << std::endl;
  }

  //DEBUG
  /**
  std::string master_debug_posedump = "master_debug_posedump.out";
  std::ofstream debug_master;
  debug_master.open(master_debug_posedump.c_str(),std::fstream::app);

  for( core::Size ii = 1; ii <= Pool_RMSD::size(); ii++ ) {
    FArray2D<double> tmp;
    std::string tag = Pool_RMSD::get_tag( ii );
    Pool_RMSD::get( ii, tmp );
    //write to file
    debug_master << tag << "  BEGIN " << std::endl;
    for( core::Size f_index_i = 1; f_index_i <= tmp.u2(); f_index_i++ ){
      for( core::Size f_index_j = 1; f_index_j <= tmp.u1(); f_index_j++ ){
	debug_master << tmp( f_index_j, f_index_i ) << " ";
      }
      debug_master << std::endl;
    }
    debug_master << tag << " END " << std::endl;
  }
  **/
  //DEBUG

  return;
}


bool MPIBPool_RMSD::is_master_node(){
  return rank_ == master_node_;
}

core::Size MPIBPool_RMSD::evaluate_and_add(
  core::pose::Pose const& pose,
  std::string& best_decoy,
  core::Real& best_rmsd,
  core::Real transition_threshold
  ){
  tr << "size of pool is " << Pool_RMSD::size() << std::endl;
  PROF_STOP( basic::MPI_POOL_SLAVE_THINKS );
  //bool use_broadcasting = true;
  assert(transition_threshold == transition_threshold_);
  if( transfer_buf_.nresidues_ == 0 ){
    transfer_buf_.nresidues_ = pose.size();
  }else{
    assert( transfer_buf_.nresidues_ == pose.size() );
  }
  if( tracer_visible_ ){
    tr.Debug << "this node is rank " << rank_ << " pool-rank is " << pool_rank_ << " master node is rank " << master_node_ << " and total size is " << npes_ << " pool-size is " << pool_npes_ << std::endl;
  }
  assert(pool_rank_ > pool_master_node_ && pool_rank_ < pool_npes_);
  //tr.Debug << "node is rank " << rank_ << " out of " << npes_ << std::endl;
  core::Size best_index;
  best_index = Pool_RMSD::evaluate( pose, best_decoy, best_rmsd );
  if( tracer_visible_ ){
    tr.Debug << "best rmsd after first evaluation is " << best_rmsd << " threashold " << transition_threshold << " and index  " << best_index << std::endl;
  }

  std::string new_tag;
  if( best_rmsd > transition_threshold ){
    //FArray2D_double coords( 3, pose.size() , 0.0 );
    PROF_START( basic::FARRAY_MANIPULATION );
    reformat( pose, new_tag );
    PROF_STOP( basic::FARRAY_MANIPULATION );
    if( tracer_visible_ ){
      tr.Debug << " slave: about to report new coordinates to master!" << std::endl;
    }

    slave_gather_new_coords();
  }else{
    if( tracer_visible_ ){
      tr.Debug << " slave: about to report, no new coordinates, to master!" << std::endl;
    }
    slave_report_no_new_coords();
  }

  int new_size = 0;
  int num_structures_to_add = 0;
  PROF_START( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );
  MPI_Bcast( &num_structures_to_add, 1, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
  PROF_STOP( basic::MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD );

  if( tracer_visible_ ){
    tr.Debug << " expecting " << num_structures_to_add << " from master node " << std::endl;
  }
  PROF_START( basic::MPI_MASTER_BCAST_WINNING_RANKS );
  MPI_Bcast( transfer_buf_.winning_ranks_, num_structures_to_add, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
  PROF_STOP( basic::MPI_MASTER_BCAST_WINNING_RANKS );
  if( tracer_visible_ ){
    tr.Debug << "received "<< num_structures_to_add << " from the master node! " << std::endl;
    for ( core::Size ii = 0; static_cast<int>(ii) < num_structures_to_add; ii++ ) {
      tr.Debug << "winning rank: " << transfer_buf_.winning_ranks_[ ii ] << "\n";
    }
    tr.Debug << std::endl;
  }
  bool i_am_a_winning_rank = false;
  for( int ii = 0; ii < num_structures_to_add; ii++ ){
    if( (int)(pool_rank_) == transfer_buf_.winning_ranks_[ ii ] ){
      if( tracer_visible_ ){
	tr.Debug << "I WON! I'm one of the winning ranks! " << std::endl;
      }
      i_am_a_winning_rank = true;
      best_rmsd = 0.0;
      assign_tag( new_tag, ( Pool_RMSD::size() + ii ) );
      tr.Debug << "I'm gonna add a new structure, this is it's tag: " << new_tag << std::endl;
      best_decoy = new_tag;
      best_decoy = new_tag;
      PROF_START( basic::WRITE_TO_FILE );
      core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();

      ss->fill_struct( pose, new_tag );
      core::io::silent::SilentFileData sfd;
      sfd.write_silent_struct( *ss, new_decoys_out_, false );
      PROF_STOP( basic::WRITE_TO_FILE );
    }
  }

  broadcast_newest_coords( num_structures_to_add ); //automatically increments pool_size appropriately
  PROF_START( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );
  MPI_Bcast( &new_size, 1, MPI_INT, pool_master_node_, MPI_COMM_POOL );
  PROF_STOP( basic::MPI_MASTER_BCAST_NEW_COMM_SIZE );

  if( tracer_visible_ ){
    tr.Debug << "new size of pool " << new_size << " current_size of pool_npes_ " << pool_npes_ << std::endl;
  }
  if( new_size != static_cast<int>(pool_npes_) ){

    PROF_START( basic::MPI_MASTER_BCAST_NEW_POOL_RANKS );
    MPI_Bcast( transfer_buf_.int_buf1_, new_size, MPI_INT, (pool_master_node_), MPI_COMM_POOL );
    PROF_STOP( basic::MPI_MASTER_BCAST_NEW_POOL_RANKS );
    if( tracer_visible_ ){
      tr.Debug << "creating new communicator from ranks: ";
      for(int ii = 0; ii < new_size; ii++){
	tr.Debug << transfer_buf_.int_buf1_[ ii ] << " ";
      }
      tr.Debug << std::endl;
    }
    create_comm( transfer_buf_.int_buf1_, new_size  );
  }
  PROF_START( basic::MPI_POOL_SLAVE_THINKS );

  if( !i_am_a_winning_rank  && num_structures_to_add > 0 ){ //update rms info and nearest cluster info
    FArray2D_double coords( 3, pose.size(), 0.0 );
    toolbox::fill_CA_coords( pose, coords );
    tr.Debug << "checking coords" << std::endl;

    tr.Debug <<std::endl;
    core::Real competing_best_rmsd= -1;
    std::string competing_best_decoy;
    if( tracer_visible_ ){
      tr.Debug << "before update, this is my info: " << best_rmsd << " " << best_decoy << " " << best_index << std::endl;
    }
    core::Size alt_index = Pool_RMSD::evaluate( coords, competing_best_decoy, competing_best_rmsd,  ( Pool_RMSD::size() - num_structures_to_add + 1 ) );
    if( tracer_visible_ ){
      tr.Debug << "after 2nd eval, this is my info: " << competing_best_rmsd << " " << competing_best_decoy << " " << alt_index << std::endl;
    }
    if( competing_best_rmsd < best_rmsd ) {
      best_rmsd = competing_best_rmsd;
      best_decoy = competing_best_decoy;
      best_index = alt_index + Pool_RMSD::size() - num_structures_to_add;
    }
  }
  if( tracer_visible_ ){
    tr.Debug << "best rmsd after evaluation is " << best_rmsd << " threashold " << transition_threshold << " num_structures_to_add " << num_structures_to_add << " pool-size " << Pool_RMSD::size() << " and index " << ( Pool_RMSD::size() - num_structures_to_add + 1 ) << std::endl;
  }

  return best_index;
}


} //mc_convergence_checks
} //moves
} //protocols
#else
#endif
