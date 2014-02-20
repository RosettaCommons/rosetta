// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIBPool_ConvergenceCheck_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIBPool_ConvergenceCheck_hh
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/MPIBPool_ConvergenceCheck.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

//Auto Headers



#ifdef USEMPI
#include <mpi.h>
#endif


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

typedef ObjexxFCL::FArray2D<double> FArray2D_double;
typedef ObjexxFCL::FArray3D<double> FArray3D_double;

struct TransferBuffer {
public:
  TransferBuffer();
  TransferBuffer( core::Size num_slave_nodes );

  void set_size( int num_slave_nodes );

  ~TransferBuffer();
  //option1:
  int* memory_offset_;
  int* size_per_coords_;
  int* int_buf1_;
  int* winning_ranks_;
  double* farray_coord_ptr_;
  FArray2D_double temp_coord_for_evaluation_;
  FArray3D_double coords_;
  core::Size size_; // num-structures reported, to broadcast
  core::Size nresidues_;
};


class MPIBPool_RMSD : public Pool_RMSD {

public:

static int master_node_;
static int pool_master_node_;



MPIBPool_RMSD( std::string const& silent_file );


void create_comm( int  ranks_to_include[], int new_size );

void update_ranks( int const active_nodes[], int new_size );

void set_discovered_out( std::string const& newout);

std::string const& get_discovered_out();

bool is_active_node();

void partition_into_coordinates(
				double const received_array[],
				int const size_per_coord[],
				core::Size num_received,
				utility::vector1<FArray2D_double>& coords
				);

void farray_to_array( core::Size index );

void farray_to_array( core::Size index, core::Size num_to_add );

void array_to_farray( core::Size index );

void array_to_farray( core::Size index, core::Size num_to_add );

void set_transition_threshold(
  core::Real threshold
);

void set_nresidues(
  core::Size nres
);

void get_nresidues();

bool workers_finished();

  /**
void add_pose_to_pool(
  FArray2D_double const& coords_to_add,
  std::string& tag
);
  **/

void add_pose_to_pool();

void finalize();

void master_go();

bool is_master_node();


core::Size evaluate_and_add(
  core::pose::Pose const& pose,
  std::string& best_decoy,
  core::Real& best_rmsd,
  core::Real transition_threshold
);

private:

void initialize();

void increment_pool_size( core::Size new_structures );

void reformat( core::pose::Pose const& pose, std::string& new_tag );

void assign_tag( std::string& new_tag, core::Size optional_id_num );


void broadcast_newest_coords( int num_to_send );
void broadcast_coords( FArray2D_double & coords );
  //void master_gather_new_coords( utility::vector1<FArray2D_double> & coords, utility::vector1<int> & winning_ranks );
void master_gather_new_coords();

void slave_gather_new_coords();
void slave_report_no_new_coords();



private:
  core::Size workers_finished_;
  utility::vector1< bool > nodes_finished_;
  core::Size pool_size_;
  core::Size new_structures_;
  int rank_;
  int pool_rank_;
  int npes_;
  int pool_npes_;
  core::Real transition_threshold_;
  std::string new_decoys_out_;
  bool tracer_visible_;
  //
#ifdef USEMPI
  static MPI_Comm MPI_COMM_POOL;
#endif

  TransferBuffer transfer_buf_;

};



} //mc_convergence_checks
} //moves
} //protocols

#endif
