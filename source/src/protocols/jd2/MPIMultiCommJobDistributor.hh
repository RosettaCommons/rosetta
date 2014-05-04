// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  intended to do JD2-based distribution of jobs where each job can run a parallel MPI protocol
/// @brief  prominent example... ReplicaExchange... run N trajectories of M Replica on an NxM+2 allocation
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_MPIMultiCommJobDistributor_hh
#define INCLUDED_protocols_jd2_MPIMultiCommJobDistributor_hh

#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/MPIMultiCommJobDistributor.fwd.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

///@details This JobDistributor is intended for JD2 controlled runs of multiple parallel simulations that use
/// multiple replicas: e.g., N trajectories of M replicas on N*M+2 processes. (2 processes for central JD2 and MPIFileBuf )
/// derived from MPIFileBufJobDistributor: two dedicated processes are used to handle JobDistribution and File-IO.
/// all other processes (higher rank ) are used for computation.
class MPIMultiCommJobDistributor : public MPIFileBufJobDistributor {
	typedef MPIFileBufJobDistributor Parent;
protected:
  ///@brief ctor is protected; singleton pattern
  MPIMultiCommJobDistributor( core::Size sub_size );

	virtual void handle_interrupt() {}

public:

	///@brief dummy for master/slave version
  virtual core::Size get_new_job_id();

	///@brief overloaded to suppress message from higher-rank replicas
	virtual	void job_succeeded(core::pose::Pose & pose, core::Real run_timei, std::string const & tag);

	///@brief overloaded to suppress message from higher-rank replicas
	virtual	void job_failed( core::pose::Pose & pose, bool );

  #ifdef USEMPI
	///@brief return current communicator ( or MPI_COMM_WORLD ).
	/// this inteface could be put into an abstract base derived directly from JobDistributor
	/// JobDistributor->MPIJobDistributor->MPIXXXJobDistributor
	/// but for now it is only used by MPIMultiCommJobDistributor...
	/// access via utility function jd2::current_mpi_comm()
	MPI_Comm const& current_mpi_comm();
	#endif

	core::Size sub_rank() {
		return sub_rank_;
	}

  friend class JobDistributorFactory; //ctor access

private:

	void setup_sub_communicators( core::Size sub_size );

	core::Size n_comm_;
#ifdef USEMPI
	utility::vector1< MPI_Comm > mpi_communicators_;
	utility::vector1< MPI_Group > mpi_groups_;
	utility::vector1< int* > mpi_ranks_;
#endif /*USEMPI*/

	int sub_rank_;
	core::Size communicator_handle_;

};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_MPIMultiCommJobDistributor_HH
