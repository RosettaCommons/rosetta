#ifdef USEMPI
#include <protocols/canonical_sampling/mc_convergence_checks/MPIPool_ConvergenceCheck.hh>
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
#include  "mpi.h"

//SilentFileStuff
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

//option stuff
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>


//profiling
#include <basic/prof.hh>
//OPT_2GRP_KEY(String, sampling, out, new_structures )

namespace protocols {
namespace canonical_sampling{
namespace mc_convergence_checks {

static thread_local basic::Tracer tr( "MPIPool_ConvergenceCheck" );

  //bool protocols::canonical_sampling::mc_convergence_checks::MPIPool_RMSD::options_registered_( false );

  //
  core::Size MPIPool_RMSD::master_node_;

using namespace ObjexxFCL;
using namespace core;
using namespace utility::io::mpi_stream;
using namespace core::io::silent;

int const POSE_TRANSFER = 2000;
int const TAG_TRANSFER = 2001;
int const MPI_UPDATE = 1000;
int const MPI_ADD_POSE_TO_POOL = 1001;
int const MPI_DATA_TRANSFER = 1003;
int const MPI_BOOL = 1006;
int const UPDATE_MISSING_POSE = 1004;
int const MPI_REPORT_SIZE = 1007;
int const MPI_FINISHED = 1008;
  //int const MASTER_NODE = master_node_;



typedef FArray2P<double> FArray2P_double;
typedef FArray2D<double> FArray2D_double;



MPIPool_RMSD::MPIPool_RMSD( std::string silent_file ):
  Pool_RMSD( silent_file ),
  trajectories_finished_( 0 ),
  new_structures_( 0 ),
  rank_( 0 ),
  npes_( 0 ),
  transition_threshold_(-1),
  new_decoys_out_( "discovered_decoys.out" )
{
  initialize();
  tr.Debug << "checking: rank " << rank_ << " has " << Pool_RMSD::size() << " structures in its pool " << std::endl;
}

  /**
void
MPIPool_RMSD::register_options(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  if( !options_registered_ ){
    tr.Debug << "registering options for sampling::out::new_structures " << std::endl;
    NEW_OPT( sampling::out::new_structures, "write structures above transition_threshold to this file", "discovered_decoys.out" );
    options_registered_ = true;
  }
}

void
MPIPool_RMSD::set_defaults_from_cmdline(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  runtime_assert( options_registered_ );
  new_decoys_out_ =  option[ sampling::out::new_structures ] ;
}
  **/


void MPIPool_RMSD::initialize(){
  PROF_START( basic::CHECK_COMM_SIZE  );
  MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank_ ) );
  MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &npes_ ) );
  PROF_STOP( basic::CHECK_COMM_SIZE  );

  //  register_options();
  //  set_defaults_from_cmdline();

  if( rank_ == master_node_ ){
    tr.Debug << "initializing master pool. setting pool_size_ " << Pool_RMSD::size() << std::endl;
    pool_size_ = Pool_RMSD::size();
    Pool_RMSD::clear();
    runtime_assert( Pool_RMSD::size() == 0 );
  }else{
    pool_size_ = Pool_RMSD::size();
  }
}

void MPIPool_RMSD::increment_pool_size( core::Real num_to_add ){
  new_structures_ += num_to_add;
}

void MPIPool_RMSD::send_update( int const receiving_rank, int message_type){
  PROF_START( basic::MPI_SEND_UPDATE );
  tr.Debug << rank_ << " is sending an update message of type " << message_type << " to " << receiving_rank << std::endl;
  MPI_Send( &message_type, 1, MPI_INT, receiving_rank, MPI_UPDATE, MPI_COMM_WORLD );
  PROF_STOP( basic::MPI_SEND_UPDATE );
}

void MPIPool_RMSD::field_message( int& sending_rank, int& message_type ){

  PROF_START( basic::MPI_SEND_UPDATE );
  MPI_Status stat;
  MPI_Recv( &message_type, 1, MPI_INT, MPI_ANY_SOURCE, MPI_UPDATE, MPI_COMM_WORLD, &stat );
  sending_rank = stat.MPI_SOURCE;
  tr.Debug << "from " << sending_rank << " received an update message_type: " << message_type << std::endl;
  PROF_STOP( basic::MPI_SEND_UPDATE );
}


  // function only used by master node
void MPIPool_RMSD::send_newest_xyz( core::Size num_to_send, int const receiving_rank ){
  core::Size current_size = new_structures_;
  runtime_assert( new_structures_ >= num_to_send );
   while( num_to_send > 0 ){
    FArray2D<double> coords;
    Pool_RMSD::get( current_size - num_to_send + 1, coords );
    std::string tag = Pool_RMSD::get_tag( current_size - num_to_send + 1);
    tr.Debug << "Sending " << tag << " from rank " << rank_ <<
      " at index " << (current_size - num_to_send + 1) <<
      " current size is " << current_size <<
      " num to send is " << num_to_send <<
      " to receiving_rank " << receiving_rank << std::endl;
    send_xyz( coords, tag, receiving_rank );
    num_to_send--;
  }
}

void MPIPool_RMSD::send_xyz(
  FArray2D<double>& xyz,
  std::string& tag,
  core::Size receiving_rank
){
  std::string coords;
  farray_to_string( xyz, coords );
  std::ostringstream double_to_string;
  double_to_string.str( coords );

  PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
  int buf[ 4 ];
  buf[ 0 ] = xyz.u1();
  buf[ 1 ] = xyz.u2();
  buf[ 2 ] = (double_to_string.str()).size();
  buf[ 3 ] = tag.size();

  MPI_Send( &buf, 4, MPI_INT, receiving_rank, POSE_TRANSFER, MPI_COMM_WORLD );
  PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );
  PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  MPI_Send(
	   const_cast<char*> (double_to_string.str().data()),
	   (double_to_string.str()).size(),
	   MPI_CHAR,
	   receiving_rank,
	   MPI_DATA_TRANSFER,
	   MPI_COMM_WORLD
  );
  tr.Debug << receiving_rank << " sending tag " << tag << std::endl;
  MPI_Send( const_cast<char*> (tag.data()), tag.size(), MPI_CHAR, receiving_rank, TAG_TRANSFER, MPI_COMM_WORLD );
  PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );
}

void MPIPool_RMSD::farray_to_string(
  FArray2D<double>& xyz,
  std::string& string
){
  PROF_START( basic::FARRAY_MANIPULATION );
  std::ostringstream double_to_string;
  for( unsigned i = 1; static_cast<int>(i) <= xyz.u1(); i++ ){
    for( unsigned j = 1; static_cast<int>(j) <= xyz.u2(); j++ ){
      double_to_string << xyz( i, j ) << " ";
    }
  }
  string = double_to_string.str();
  PROF_STOP( basic::FARRAY_MANIPULATION );
}

void MPIPool_RMSD::receive_newest_xyz( core::Size num_to_get, int const sending_rank ){
  while( num_to_get > 0 ){
    FArray2D<double> coords;
    std::string tag;
    receive_xyz( coords, tag, sending_rank );
    num_to_get--;
    Pool_RMSD::add( coords, coords.u2(), tag );
    tr.Debug << "added pose (receive_newest_xyz) now has " << Pool_RMSD::size() << std::endl;
  }
}

void MPIPool_RMSD::receive_xyz(
  FArray2D<double>& xyz,
  std::string& tag,
  core::Size sending_rank
){

  PROF_START( basic::MPI_SLAVE_REPORT_SIZES );
  int buf[ 4 ];
  MPI_Status stat;
  MPI_Recv( &buf, 4, MPI_INT, sending_rank, POSE_TRANSFER, MPI_COMM_WORLD, &stat );
  core::Size xyz_u1;
  core::Size xyz_u2;
  core::Size received_string_size;
  core::Size tag_size;
  xyz_u1 = buf[ 0 ];
  xyz_u2 = buf[ 1 ];
  received_string_size = buf[ 2 ];
  tag_size = buf[ 3 ];
  PROF_STOP( basic::MPI_SLAVE_REPORT_SIZES );
  //core::Size index = 0;
  xyz.dimension( xyz_u1, xyz_u2, 0.0 );
  std::string data;
  PROF_START( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  char *matrix = new char[ received_string_size + 1 ];

  MPI_Recv( matrix, received_string_size , MPI_CHAR, sending_rank, MPI_DATA_TRANSFER, MPI_COMM_WORLD, &stat );
  data.assign( matrix, received_string_size );
  string_to_farray( xyz, data, xyz_u1, xyz_u2 );
  char *cbuf = new char[tag_size+1];
  MPI_Recv( cbuf, tag_size, MPI_CHAR, sending_rank, TAG_TRANSFER, MPI_COMM_WORLD, &stat );
  tag.assign( cbuf, tag_size );
  PROF_STOP( basic::MPI_SLAVE_REPORT_NEW_COORDS );
  //tr.Debug << sending_rank << " RECEIVED TAG " << tag << std::endl;
}

void MPIPool_RMSD::string_to_farray( FArray2D<double>& xyz, std::string& string, int xyz_u1, int xyz_u2 ){
  std::istringstream string_to_double;
  string_to_double.str(string);
  PROF_START( basic::FARRAY_MANIPULATION );
  double element;
  for(core::Size i = 1; static_cast<int>(i) <= xyz_u1; i++ ){
    for(core::Size j = 1; static_cast<int>(j) <= xyz_u2; j++ ){
      string_to_double >> element;
      xyz( i, j ) = element;
    }
  }
  PROF_STOP( basic::FARRAY_MANIPULATION );
}

void MPIPool_RMSD::set_transition_threshold( core::Real threshold ){
  transition_threshold_ = threshold;
}

void MPIPool_RMSD::set_discovered_out( std::string new_out ){
  new_decoys_out_ = new_out;
}

std::string MPIPool_RMSD::get_discovered_out(){
  return new_decoys_out_;
}

core::Size MPIPool_RMSD::get_pool_diff(core::Size target_rank){
  int target_rank_pool_size;
  int my_pool_size = pool_size_ + new_structures_;
  tr.Debug << "(get_pool_diff) current size is " << my_pool_size << std::endl;
  int num_diff;
  MPI_Status status;
  if(rank_ == master_node_){
    PROF_START( basic::MPI_SYNC_POOL_DIFF );
    MPI_Recv( &target_rank_pool_size, 1, MPI_INT, target_rank, MPI_REPORT_SIZE, MPI_COMM_WORLD, &status );
    num_diff = my_pool_size - target_rank_pool_size;
    MPI_Send( &num_diff, 1, MPI_INT, target_rank, MPI_REPORT_SIZE, MPI_COMM_WORLD );
    tr.Debug << "MASTER NODE HAS " << my_pool_size << " AND SLAVE HAS " << target_rank_pool_size << " SO NEED TO UPDATE " << num_diff << " LAST POSES " << std::endl;
    PROF_STOP( basic::MPI_SYNC_POOL_DIFF );
  }else{
    PROF_START( basic::MPI_SYNC_POOL_DIFF );
    MPI_Send( &my_pool_size, 1, MPI_INT, master_node_, MPI_REPORT_SIZE, MPI_COMM_WORLD );
    MPI_Recv( &num_diff, 1, MPI_INT, master_node_, MPI_REPORT_SIZE, MPI_COMM_WORLD, &status );
    PROF_STOP( basic::MPI_SYNC_POOL_DIFF );
  }

  return num_diff;
}

void MPIPool_RMSD::send_accepted(bool truefalse, core::Size rank){
  int bool_to_int = 0;
  if( truefalse ){
    bool_to_int = 1;
  }
  tr.Debug << "(send_accepted)sending to rank " << rank << " " << truefalse << std::endl;
  PROF_START( basic::MPI_SEND_ACCEPTED );
  MPI_Send( &bool_to_int, 1, MPI_INT, rank, MPI_BOOL, MPI_COMM_WORLD );
  PROF_STOP( basic::MPI_SEND_ACCEPTED );
}

bool MPIPool_RMSD::receive_is_accepted(core::Size rank){
  int bool_to_int;
  MPI_Status stat;
  PROF_START( basic::MPI_SEND_ACCEPTED );
  MPI_Recv( &bool_to_int, 1, MPI_INT, rank, MPI_BOOL, MPI_COMM_WORLD, &stat );
  PROF_STOP( basic::MPI_SEND_ACCEPTED );
  if( bool_to_int == 1){
    return true;
  }else{
    return false;
  }
}


bool MPIPool_RMSD::trajectories_finished(){
  if(trajectories_finished_ < ( npes_ - master_node_ - 1 ) ){
    tr.Debug << "num trajectories finished: " <<
      trajectories_finished_ << " needed: " <<
      ( npes_ - master_node_ - 1 ) << std::endl;
    return false;
  }else{
    tr.Debug << "FINISHED! num trajectories finished: " <<
      trajectories_finished_ << " needed: " <<
      ( npes_ - master_node_ - 1 ) << std::endl;
    return true;
  }
}

void MPIPool_RMSD::finalize(){
  if( rank_ != master_node_ ){
    tr.Debug << "sending finalized message to master" << std::endl;
    send_update( master_node_, MPI_FINISHED );
  }
}


void MPIPool_RMSD::master_go(){
  tr.Debug << "in master go master rank is " << rank_ << std::endl;

  PROF_START( basic::MPICANONICALSAMPLING );
  while( !trajectories_finished() ){

    int slave_rank = -1;
    int message_type;
    field_message( slave_rank /**receive update from any node*/, message_type );

    if( message_type == UPDATE_MISSING_POSE ){
      core::Size num_to_update = get_pool_diff( slave_rank );
      send_newest_xyz( num_to_update, slave_rank );
    }
    else if(message_type == MPI_ADD_POSE_TO_POOL){
      //tr.Debug << "received message to add pose to pool" << std::endl;
      FArray2D<double> new_coords;
      std::string new_tag;

      receive_xyz( new_coords, new_tag, slave_rank );
      core::Size num_recent_updates = get_pool_diff( slave_rank );
      send_newest_xyz( num_recent_updates, slave_rank );
      core::Real best_rmsd = -1; std::string best_decoy;
      bool is_accepted = false;
      if( num_recent_updates == 0 ){ //no additions during evaluation
	Pool_RMSD::add( new_coords, new_coords.u2(), new_tag );
	increment_pool_size(1);
	// DEBUG OUTPUT
	/**
	std::ofstream debug_cl_cnters;
	std::ostringstream q;
	q << rank_;
	debug_cl_cnters.open((q.str() + ".debug_cl_centers.txt").c_str(),std::fstream::app);
	  debug_cl_cnters << new_tag << std::endl;
	**/
	//DEBUG OUTPUT

	  tr.Debug << "MASTER: coords just added. pool size is now " << (pool_size_ + new_structures_)
		   << " pool_size_ is " << pool_size_
		   << " new_structures_ " << new_structures_
		   << " real size is: " << Pool_RMSD::size() << std::endl;
	is_accepted = true;
      }else{
	runtime_assert( new_structures_ >= num_recent_updates );
	tr.Debug << "starting evaluation from index: " << (new_structures_ - num_recent_updates +1) << " of " << new_structures_ << std::endl;
	Pool_RMSD::evaluate( new_coords, best_decoy, best_rmsd, (new_structures_ - num_recent_updates +1) );
	runtime_assert( transition_threshold_ != -1 );
	if( best_rmsd > transition_threshold_ ){
	  tr.Debug << "best_rmsd " << best_rmsd <<
	  	    " greater than transition threshold " << transition_threshold_ << std::endl;
	  is_accepted = true;
	  Pool_RMSD::add( new_coords, new_coords.u2(), new_tag );
	  increment_pool_size(1);
	  tr.Debug << " coords just added. pool size is now " << (pool_size_ + new_structures_) << " real size is " << Pool_RMSD::size() << std::endl;
	  // DEBUG OUTPUT
	  /**
	  std::ofstream debug_cl_cnters;
	  std::ostringstream q;
	  q << rank_;
	  debug_cl_cnters.open((q.str() + ".debug_cl_centers.txt").c_str(),std::fstream::app);
	  debug_cl_cnters << new_tag << std::endl;
	  **/
	  //DEBUG OUTPUT
	}
      }
      send_accepted( is_accepted, slave_rank );
   }else if( message_type == MPI_FINISHED ){
      trajectories_finished_++;
    }
  }
  PROF_STOP( basic::MPICANONICALSAMPLING );
  tr.Debug << "master node finished " << std::endl;
  return;
}

bool MPIPool_RMSD::is_master_node(){
  tr.Debug << "testing if is master node:  rank is " << rank_ << " and master node rank is " << master_node_ << std::endl;
  return rank_ == master_node_;
}

core::Size MPIPool_RMSD::evaluate_and_add(
  core::pose::Pose const& pose,
  std::string& best_decoy,
  core::Real& best_rmsd,
  core::Real transition_threshold
  ){

  runtime_assert(transition_threshold == transition_threshold_);

  tr.Debug << "node is rank " << rank_ << " out of " << npes_ << std::endl;
  tr.Debug << " using MPIPool_RMSD::evaluate_and_add" << std::endl;
  core::Size best_index = -1;

  PROF_START( basic::MPICANONICALSAMPLING );

  send_update( master_node_, UPDATE_MISSING_POSE );
  core::Size num_to_update = get_pool_diff( rank_ );
  receive_newest_xyz( num_to_update, master_node_ );
  increment_pool_size( num_to_update );
  tr.Debug << " slave pool just received " << num_to_update << " pool_size is now " << (pool_size_ + new_structures_ ) << " real size: " << Pool_RMSD::size() << std::endl;
  //runtime_assert(get_pool_diff( rank_ ) == 0);

  PROF_STOP( basic::MPICANONICALSAMPLING );

  best_index = Pool_RMSD::evaluate( pose, best_decoy, best_rmsd );
  tr.Debug << "best rmsd after evaluation is " << best_rmsd << " threadhol " << transition_threshold << std::endl;
  if(best_rmsd > transition_threshold){
    tr.Debug << "best_rmsd is " << best_rmsd << " which is greater than transition_threshold, adding pose to pool " << std::endl;

    PROF_START( basic::MPICANONICALSAMPLING );

    send_update( master_node_, MPI_ADD_POSE_TO_POOL );
    FArray2D_double coords( 3, pose.total_residue() , 0.0 );
    PROF_START( basic::FARRAY_MANIPULATION );
    protocols::toolbox::fill_CA_coords( pose, pose.total_residue(), coords );
    PROF_STOP( basic::FARRAY_MANIPULATION );
    //assign new tag based on olli's scheme
    std::string jobname = protocols::jd2::current_output_name();
    std::string new_cluster_tag = "new."+lead_zero_string_of( Pool_RMSD::size(), 8 )+".0"+"_"+jobname;

    send_xyz( coords, new_cluster_tag, master_node_ );
    core::Size num_recent_updates = get_pool_diff( rank_ );
    receive_newest_xyz( num_recent_updates, master_node_ );
    increment_pool_size( num_recent_updates );
    tr.Debug << " slave pool just received " << num_recent_updates << " pool_size is now " << ( pool_size_ + new_structures_ ) << " real size: " << Pool_RMSD::size() << std::endl;
    //runtime_assert(get_pool_diff( rank_ ) == 0);
    bool is_accepted = receive_is_accepted( master_node_ );

    PROF_STOP( basic::MPICANONICALSAMPLING );

    if( is_accepted ){
      tr.Debug << "master accepted new pose(evaluate). adding new pose to pool, "
	       << Pool_RMSD::size() << std::endl;

      Pool_RMSD::add( coords, coords.u2(), new_cluster_tag );

      PROF_START( basic::WRITE_TO_FILE );
      core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
      ss->fill_struct( pose, new_cluster_tag );
      core::io::silent::SilentFileData sfd;
      sfd.write_silent_struct( *ss, new_decoys_out_, false );
      PROF_STOP( basic::WRITE_TO_FILE );

      increment_pool_size(1);
      tr.Debug << " coords just added size is now " << ( pool_size_ + new_structures_ ) << " and real size is " << Pool_RMSD::size() << std::endl;
    }
  }
  return best_index;
}


} //mc_convergence_checks
} //moves
} //protocols
#else
#endif
