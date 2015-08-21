// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/MPI_Master.hh
/// @brief  The MPI_Master role in elscripts, handles trajectories, generating workunits, processing of results
/// Has additional functions related to boost mpi world and MPI_EndPoint
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_MPI_Master_hh
#define INCLUDED_protocols_elscripts_MPI_Master_hh
#if defined (USEBOOSTMPI) && defined (USELUA)
// this is useless without mpi and lua
#include <protocols/elscripts/MPI_Master.fwd.hh>
#include <protocols/elscripts/Master.hh>
#include <protocols/wum2/MPI_EndPoint.hh>

namespace protocols {
namespace elscripts {

void lregister_MPI_Master( lua_State * lstate );

class MPI_Master : public Master {
public:
	// default memory limit is 1GB
	// default reserved mem size is 100MB as recommended by fpd
	MPI_Master( boost::mpi::communicator world, std::vector<int> slaves, int num_trajectories = 1, boost::uint64_t mem_limit=2147483648, boost::uint64_t reserved_mem=104857600, boost::uint64_t reserved_mem_multiplier=5 );
	~MPI_Master(){}
	void go();

private:
	int inputter_rank() {
		// handles if there is or is not pool logic, needed for inputter offset
		// doesnt do anything now
		return world_.rank() + 1;
	}

	//void request_pool_structures( std::vector< int > needs_replace );
	//void request_pool_structure( int trajectory_idx );

private:
	boost::mpi::communicator world_;

	std::vector< int > slaves_;
	boost::posix_time::ptime last_status_sweep_time_;
};

} //elscripts
} //protocols
#endif
#endif
