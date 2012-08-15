// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/Master.hh
/// @brief  The Master role in elscripts, handles trajectories, generating workunits, processing of results
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_Master_hh
#define INCLUDED_protocols_elscripts_Master_hh
#if defined (USEBOOSTMPI) && defined (USELUA)
// this is useless without mpi and lua
#include <protocols/elscripts/Master.fwd.hh>
#include <protocols/wum2/MPI_EndPoint.hh>
#include <protocols/elscripts/BaseRole.hh>

namespace protocols {
namespace elscripts {

void lregister_Master( lua_State * lstate );

class Master : public BaseRole {
  public:
    // default memory limit is 1GB
    // default reserved mem size is 100MB as recommended by fpd
    Master( boost::mpi::communicator world, std::vector<int> slaves, int num_trajectories = 10, boost::uint64_t mem_limit=2147483648, boost::uint64_t reserved_mem=104857600, boost::uint64_t reserved_mem_multiplier=5 );
    ~Master(){}
    void go();

    boost::uint64_t available_mem() {
      boost::uint64_t buff_mem = 
        pool_comm_->current_mem() +
        slave_comm_->current_mem() +
				reserved_mem_ * reserved_mem_multiplier_ +
        trajectories_mem();
      if( buff_mem >= mem_limit_ ) {
        return 0;
      } else {
        return mem_limit_ - buff_mem;
      }
    }

		void make_wu( std::string const & wuname, int traj_idx, core::pose::Pose * p);
		void make_wu_until_limit( std::string const & wuname, int num);
		void end_traj( int traj_idx );

		void interpreter();

  private:
    boost::uint64_t current_mem() {
      return slave_comm_->current_mem() +
        pool_comm_->current_mem() +
        trajectories_mem(); // slaves_ memory trivial
    }

		int inputter_rank() {
			// handles if there is or is not pool logic, needed for inputter offset
			// doesnt do anything now
			return world_.rank() + 1;
		}

    // use Inputters to fill trajectory pipe up to limit
    // if inputters have no more decoys, then ask pool for structures
    void fill_trajectories();

    //void request_pool_structures( std::vector< int > needs_replace );
    //void request_pool_structure( int trajectory_idx );
		
    void update_trajectories_mem();
    boost::uint64_t trajectories_mem() { return trajectories_mem_; }

  private:
		boost::mpi::communicator world_;
    protocols::wum2::EndPointSP slave_comm_;
    protocols::wum2::EndPointSP pool_comm_;

    int num_trajectories_;
    std::vector< int > slaves_;
    core::io::serialization::PipeSP trajectories_;
    boost::uint64_t trajectories_mem_; 

		int num_trajectories_finished_;
};

} //elscripts
} //protocols
#endif
#endif
