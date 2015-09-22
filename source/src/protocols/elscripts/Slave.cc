// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/Slave.cc
/// @brief  the slave role of elscripts
/// @author Ken Jung

#ifdef USELUA
#include <protocols/elscripts/Slave.hh>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <protocols/wum2/WorkUnit.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace elscripts {

void lregister_Slave( lua_State * lstate ) {
	lregister_BaseRole( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<Slave, BaseRole>("Slave")
		]
	];
}

static THREAD_LOCAL basic::Tracer TR( "protocols.elscripts.Slave" );

Slave::Slave( int master, boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier) :
	master_(master),
	BaseRole( mem_limit, reserved_mem, reserved_mem_multiplier) {

		// endpoint uses this function to get role-wide memory usage
    boost::function< boost::uint64_t ()> ref_available_mem = boost::bind( &protocols::elscripts::Slave::available_mem, this );

    master_comm_ = protocols::wum2::EndPointSP( new protocols::wum2::EndPoint( ref_available_mem ) );

		lua_init();
		lregister_Slave(lstate_);
		luabind::globals(lstate_)["slave"] = this;
		luabind::globals(lstate_)["rank"] = 0;

		register_calculators();
		instantiate_tasks();
		instantiate_scorefxns();
		instantiate_filters();
		instantiate_movers();
		instantiate_output();
		instantiate_workunits();
}

void Slave::go(){
	using namespace utility::lua;

    // run the first wu in the queue
	if( ! master_comm_->inq().empty() ) {
    protocols::wum2::WorkUnitSP wu = master_comm_->inq().pop_front(); 
		// try dynamic casts to WorkUnit_elscriptsState
		if( wu ) {
			protocols::wum2::WorkUnit_ElScriptsSP castattempt = boost::dynamic_pointer_cast<protocols::wum2::WorkUnit_ElScripts> (wu);
			if( castattempt ) {
				std::string wuname = castattempt->name();
				// create temporary environment that will be thrown away after run()
				std::string action = "DELIM(\n"
					"tmp_run_env = {}\n"
					"setmetatable(tmp_run_env, {__index = _G })\n"
					")DELIM";
				int err = luaL_dostring ( lstate_, action.c_str() );
				if( err == 1) {
					TR << "Creating tmp namespace for run_on_slave() on slave failed. Error is:" << std::endl;
					TR << lua_tostring(lstate_, -1) << std::endl;
					std::exit(9);
				}

				// exporting useful stuff to the tmp env
				luabind::globals(lstate_)["tmp_run_env"]["pipemap"] = castattempt->pipemap().lock();
				luabind::globals(lstate_)["tmp_run_env"]["state"] = castattempt->state().lock();
				luabind::globals(lstate_)["tmp_run_env"]["traj_idx"] = castattempt->trajectory_idx();

				//calling run fxn
				action = "(\n"
					"els_setenv(tmp_run_env)\n"
					"els.workunits.)"+wuname+"DELIM(.run_on_slave()\n"
					"tmp_run_env = {} -- delete tmp_env after calling run()\n"
					")DELIM";
				err = luaL_dostring ( lstate_, action.c_str() );
				if( err == 1) {
					TR << "Calling lua function for workunit " << wuname << " run_on_slave fxn failed. Error is:" << std::endl;
					TR << lua_tostring(lstate_, -1) << std::endl;
					std::exit(9);
				}
				lua_gc(lstate_, LUA_GCCOLLECT, 0);
				// shallow copy! works fine for my purposes
				protocols::wum2::WorkUnitSP result_wu( new protocols::wum2::WorkUnit_ElScripts( *castattempt ) );
				master_comm_->outq().push_back( result_wu );
			} else {
				wu->run();
			}
		}
	}
}


} //elscripts
} //protocols
#endif
