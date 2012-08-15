// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/SingleNode.hh
/// @brief  the singlenode role of elscripts
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_SingleNode_hh
#define INCLUDED_protocols_elscripts_SingleNode_hh
#ifdef USELUA
// this is useless without lua
#include <lua.hpp>
#include <boost/cstdint.hpp>

#include <protocols/elscripts/SingleNode.fwd.hh>
#include <protocols/elscripts/Slave.fwd.hh>
#include <protocols/elscripts/Master.fwd.hh>

namespace protocols {
namespace elscripts {

void lregister_SingleNode( lua_State * lstate );

class SingleNode {
  public:
    // default memory limit is 2GB
    // default reserved mem size is 100MB as recommended by fpd
    SingleNode( boost::uint64_t mem_limit=2147483648, boost::uint64_t reserved_mem=104857600, boost::uint64_t reserved_mem_multiplier=10 );
    ~SingleNode(){}
    void go();

  private:
		MasterSP master_;
		SlaveSP slave_;
};

} //elscripts
} //protocols
#endif
#endif
