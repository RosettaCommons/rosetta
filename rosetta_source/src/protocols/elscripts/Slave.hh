// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/Slave.hh
/// @brief  the slave role of elscripts
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_Slave_hh
#define INCLUDED_protocols_elscripts_Slave_hh
#if defined (USEBOOSTMPI) && defined (USELUA)
// this is useless without mpi and lua
#include <protocols/elscripts/Slave.fwd.hh>
#include <protocols/wum2/EndPoint.hh>
#include <protocols/elscripts/BaseRole.hh>

namespace protocols {
namespace elscripts {

void lregister_Slave( lua_State * lstate );

class Slave : public BaseRole {
  public:
    // default memory limit is 2GB
    // default reserved mem size is 100MB as recommended by fpd
    Slave( boost::mpi::communicator world, int master, boost::uint64_t mem_limit=2147483648, boost::uint64_t reserved_mem=104857600, boost::uint64_t reserved_mem_multiplier=10 );
    ~Slave(){}
    void go();

    boost::uint64_t available_mem() {
      boost::uint64_t buff_mem =
        master_comm_->buffered_current_mem() +
        mover_cache_mem();
      if( buff_mem >= mem_limit_ ) {
        return 0;
      } else {
        return mem_limit_ - buff_mem;
      }
    }

  private: 
    boost::uint64_t current_mem() {
      return master_comm_->current_mem() + mover_cache_mem();
    }

  private:
	  boost::mpi::communicator world_;
    protocols::wum2::EndPointSP master_comm_;

    int master_;
};

} //elscripts
} //protocols
#endif
#endif
