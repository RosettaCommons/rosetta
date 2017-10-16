// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIPool_ConvergenceCheck_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIPool_ConvergenceCheck_hh

#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/MPIPool_ConvergenceCheck.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {


#ifdef USEMPI
class MPIPool_RMSD : public Pool_RMSD {

public:

static core::Size master_node_;

MPIPool_RMSD( std::string silent_file );

  /**
void register_options();

void set_defaults_from_cmdline();
  **/

void send_update(
  int const receiving_rank,
  int const message_type
);

void field_message(
  int& sending_rank,
  int& message_type
);

void farray_to_string(
  ObjexxFCL::FArray2D<double>& xyz,
  std::string& string
);

void string_to_farray(
  ObjexxFCL::FArray2D<double>& xyz,
  std::string& string,
  int xyz_u1,
  int xyz_u2
);

void set_transition_threshold(
  core::Real threshold
);

void set_discovered_out( std::string new_out );

std::string get_discovered_out();

  /**
void broadcast_new_coords(
  ObjexxFCL::FArray2D<double>& xyz,
  std::string& tag
);
  **/

bool trajectories_finished();

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
  core::Size trajectories_finished_;
  core::Size pool_size_;
  core::Size new_structures_;
  core::Size rank_;
  core::Size npes_;
  core::Real transition_threshold_;
  std::string new_decoys_out_;

  //  static bool options_registered_;


  void initialize();

  void increment_pool_size( core::Real new_structures );

  void send_newest_xyz(
    core::Size num_to_get,
    int const receiving_rank
  );

  void send_xyz(
    ObjexxFCL::FArray2D<double>& xyz,
    std::string& tag,
    core::Size rank
  );

  core::Size get_pool_diff(core::Size rank);

  void receive_newest_xyz(
    core::Size num_to_get,
    int const sending_rank
  );

  void receive_xyz(
    ObjexxFCL::FArray2D<double>& xyz,
    std::string& tag,
    core::Size rank
  );

  void send_accepted(bool truefalse, core::Size rank);
  bool receive_is_accepted(core::Size rank);

  //int pool_status();

};
#endif

} //mc_convergence_checks
} //moves
} //protocols

#endif
