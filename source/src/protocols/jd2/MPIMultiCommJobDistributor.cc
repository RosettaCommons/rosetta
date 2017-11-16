// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MPIMultiCommJobDistributor.cc
/// @brief  implementation of MPIMultiCommJobDistributor
/// @author Oliver Lange olange@u.washington.edu
/// @detail freely based on the MPIWorkPoolJobDistributor from Doug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/MPIMultiCommJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/jd2/MpiFileBuffer.hh>
#include <utility/io/ozstream.hh> //to toggle MPI rerouting

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// C++ headers
#include <string>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.jd2.MPIMultiCommJobDistributor" );

namespace protocols {
namespace jd2 {

using namespace core;


using namespace basic::options;
using namespace basic::options::OptionKeys;

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIMultiCommJobDistributor::MPIMultiCommJobDistributor( core::Size sub_size ) {
	setup_sub_communicators( sub_size );
}

void MPIMultiCommJobDistributor::setup_sub_communicators( Size sub_size ) {
	n_comm_ = ( n_rank()-min_client_rank() ) / sub_size;
	tr.Debug << " can allocate " << n_comm_ << " communication groups " << std::endl;
	tr.Debug << " n_rank: " << n_rank() << " sub_size: " << sub_size << std::endl;
	set_n_worker( n_comm_ );
#ifdef USEMPI
	if ( n_comm_ < 1 ) {
		tr.Error << "requested sub-communicators of size " << sub_size << " but only " << n_rank() << " processes are available " << std::endl;
		utility_exit_with_message( "cannot run with requested size of sub-partition" );
	}

	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group );
	communicator_handle_ = 0;
	Size i_rank = min_client_rank();
	mpi_groups_.resize( n_comm_, MPI_GROUP_NULL );
	mpi_communicators_.resize( n_comm_, MPI_COMM_NULL );
	for ( Size i_comm = 1; i_comm <= n_comm_; ++i_comm ) {
		mpi_ranks_.push_back( new int[ sub_size ] ); //delete never called because singleton class
		for ( Size i = 0; i < sub_size; ++i, ++i_rank ) {
			mpi_ranks_.back()[ i ] = i_rank;
			if ( i_rank == rank() ) {
				communicator_handle_ = i_comm;
			}
		}
		MPI_Group_incl( world_group, sub_size, mpi_ranks_.back(), &(mpi_groups_[ i_comm ]) );
		// MPI_Comm_create
		/// Note that the call is to be executed by all processes in comm,
		/// even if they do not belong to the new group. This call applies only to intra-communicators.
		MPI_Comm_create( MPI_COMM_WORLD, mpi_groups_[ i_comm ], &(mpi_communicators_[ i_comm ]) );
	}

	runtime_assert( rank() < min_client_rank() || communicator_handle_ == 0 || communicator_handle_ <= mpi_communicators_.size() );
	if ( rank() >= min_client_rank() && communicator_handle_ ) {
		MPI_Comm_rank( mpi_communicators_[ communicator_handle_ ], &sub_rank_ );
	} else {
		sub_rank_ = -1;
	}
#endif //USEMPI
}

/// @details master-job distributes job-ids as usual. sub-masters obtain a new job from master
/// all processes within a single communication context work on the same job-id. ( Bcast )
core::Size
MPIMultiCommJobDistributor::get_new_job_id() {
	if ( rank() < min_client_rank() ) {
		return Parent::get_new_job_id();
	} else {
#ifdef USEMPI
		int new_job_id( -1 );
		if ( sub_rank_ == 0 ) {
			new_job_id = Parent::get_new_job_id(); 			//this sets batch_id()
		}
		if ( sub_rank_ >= 0 ) {
			//communicate new job and batch ids to group-members...
			runtime_assert( communicator_handle_ && communicator_handle_ <= mpi_communicators_.size() );
			int mpi_buf[ 2 ];
			mpi_buf[ 0 ] = new_job_id;
			mpi_buf[ 1 ] = current_batch_id();
			MPI_Bcast( mpi_buf, 2, MPI_INT, 0, mpi_communicators_[ communicator_handle_ ] );
			new_job_id = mpi_buf[ 0 ];
			if ( sub_rank_ > 0 ) set_batch_id( mpi_buf[ 1 ] );
			runtime_assert( new_job_id >= 0 );
			return new_job_id;
		}
#endif
	} //e.g. overhang processes that didn't fit into any of the sub-groups  -- send spin-down
	return 0;
}

/// @brief dummy for master/slave version
void
MPIMultiCommJobDistributor::job_succeeded(core::pose::Pose &pose, core::Real run_time, std::string const & tag) {
	if ( sub_rank() <= 0 ) {
		Parent::job_succeeded( pose, run_time, tag);
	}
}

/// @brief dummy for master/slave version
void
MPIMultiCommJobDistributor::job_failed(core::pose::Pose &pose, bool retry ) {
	if ( sub_rank() <= 0 ) {
		Parent::job_failed( pose, retry);
	}
}

#ifdef USEMPI
MPI_Comm const& MPIMultiCommJobDistributor::current_mpi_comm() {
	runtime_assert( communicator_handle_ );
	runtime_assert( communicator_handle_ <= mpi_communicators_.size() );
	return mpi_communicators_[ communicator_handle_ ];
}
#endif

}//jd2
}//protocols
