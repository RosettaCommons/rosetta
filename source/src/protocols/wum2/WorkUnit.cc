// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum2/WorkUnit.cc
/// @brief
/// @author Ken Jung

#include <protocols/wum2/WorkUnit.hh>
#include <basic/Tracer.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#ifndef WIN_PYROSETTA
#include <windows.h>
#endif
#endif

namespace protocols {
namespace wum2 {


static thread_local basic::Tracer TR( "protocols.wum2.WorkUnit" );

#ifdef USELUA
void lregister_WorkUnit( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("wum2")
		[
			luabind::class_<WorkUnit>("WorkUnit")
		]
	];
}
void lregister_WorkUnit_Wait( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("wum2")
		[
			luabind::class_<WorkUnit_Wait>("WorkUnit_Wait")
		]
	];
}
void lregister_WorkUnit_ElScripts( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("wum2")
		[
			luabind::class_<WorkUnit_ElScripts>("WorkUnit_ElScripts")
		]
	];
}
#endif

// ---------- WorkUnit --------------
WorkUnit::WorkUnit
( core::Size master,
	core::Size trajectory_idx ) :

	master_(master),
	trajectory_idx_(trajectory_idx)
{
	unixtime_creation_ = time(NULL);
	unixtime_start_ = 0;
	unixtime_stop_ = 0;
	create_unique_id();
	prioritize_ = false;
}

void
WorkUnit::print( std::ostream & out, bool /*verbose*/ ) const {
	out << "WU_id:         " << id_    << std::endl;
	out << "WU_time_create:" << unixtime_creation_<< std::endl;
	out << "WU_time_start: " << unixtime_start_<< std::endl;
	out << "WU_time_stop:  " << unixtime_stop_<< std::endl;
}


void WorkUnit::set_run_start(){
	unixtime_start_ = time(NULL);
}

void WorkUnit::set_run_stop(){
	unixtime_stop_ = time(NULL);
}

core::Size WorkUnit::get_run_time(){
	return unixtime_stop_ - unixtime_start_;
}

void WorkUnit::create_unique_id() {
	id_ = 0;
}

// ---------- WorkUnit_Wait --------------

WorkUnit_Wait::WorkUnit_Wait
( core::Size master,
	core::Size trajectory_idx,
	long wait_time ) :
	WorkUnit( master, trajectory_idx ),
	wait_time_(wait_time)
{}

void WorkUnit_Wait::run(){
	//TR << "Waiting for " << header.extra_data_1_ << std::endl;
#ifdef _WIN32
#ifndef WIN_PYROSETTA
	Sleep( wait_time_  * 1000 );
#endif
#else
	sleep( wait_time_ );
#endif
}


// ---------- WorkUnit_ElScripts --------------

WorkUnit_ElScripts::WorkUnit_ElScripts (
	core::Size master,
	core::Size trajectory_idx,
	core::io::serialization::PipeMapSP p,
	protocols::moves::SerializableStateSP state,
	std::string name
) :
	WorkUnit( master, trajectory_idx ),
	name_(name),
	pipemap_(p),
	state_(state)
{}

void WorkUnit_ElScripts::run(){
}


} //wum2
} //protocols

