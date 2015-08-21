// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/MPI_Slave.hh
/// @brief  the slave role of elscripts
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_MPI_Slave_hh
#define INCLUDED_protocols_elscripts_MPI_Slave_hh
#if defined (USEBOOSTMPI) && defined (USELUA)
// this is useless without mpi and lua
#include <protocols/elscripts/MPI_Slave.fwd.hh>
#include <protocols/wum2/MPI_EndPoint.hh>
#include <protocols/elscripts/Slave.hh>

namespace protocols {
namespace elscripts {

void lregister_MPI_Slave( lua_State * lstate );

class MPI_Slave : public Slave {
public:
	// default memory limit is 2GB
	// default reserved mem size is 100MB as recommended by fpd
	MPI_Slave( boost::mpi::communicator world, int master, boost::uint64_t mem_limit=2147483648, boost::uint64_t reserved_mem=104857600, boost::uint64_t reserved_mem_multiplier=10 );
	~MPI_Slave(){}
	void go();

private:
	boost::mpi::communicator world_;
};

} //elscripts
} //protocols
#endif
#endif
