// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/SingleNode.cc
/// @brief  the slave role of elscripts
/// @author Ken Jung

#ifdef USELUA
// this is useless without lua
#include <protocols/elscripts/SingleNode.hh>

#include <protocols/elscripts/Slave.hh>
#include <protocols/elscripts/Master.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace elscripts {

void lregister_SingleNode( lua_State * lstate ) {
	lregister_Slave( lstate );
	lregister_Master( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<SingleNode>("SingleNode")
		]
	];
}

static thread_local basic::Tracer TR( "protocols.elscripts.SingleNode" );

SingleNode::SingleNode( boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier) {
	// here they share the same memory limit, which is wrong
	// they should be passed functors to SingleNode's total memory, including both roles
	// or just give them half the total mem each
	master_ = MasterSP( new Master( 1, mem_limit, reserved_mem, reserved_mem_multiplier ) );
	slave_ = SlaveSP( new Slave( 1, mem_limit, reserved_mem, reserved_mem_multiplier ) );
}

void SingleNode::go(){
  while( 1 ) {
		master_->go();
		slave_->inq().push_back(master_->outq().pop_all());
		slave_->go();
		master_->inq().push_back(slave_->outq().pop_all());
  }
}


} //elscripts
} //protocols
#endif
