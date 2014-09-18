// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/MPI_Slave.cc
/// @brief  the slave role of elscripts
/// @author Ken Jung

#if defined (USEBOOSTMPI) && defined (USELUA)
// this is useless without mpi
#include <protocols/elscripts/MPI_Slave.hh>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <protocols/wum2/WorkUnit.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace elscripts {

void lregister_MPI_Slave( lua_State * lstate ) {
	lregister_Slave( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<MPI_Slave, Slave>("MPI_Slave")
		]
	];
}

static thread_local basic::Tracer TR( "protocols.elscripts.MPI_Slave" );

MPI_Slave::MPI_Slave( boost::mpi::communicator world, int master, boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier) :
	world_(world),
	Slave( master, mem_limit, reserved_mem, reserved_mem_multiplier) {

		// endpoint uses this function to get role-wide memory usage
    boost::function< boost::uint64_t ()> ref_available_mem = boost::bind( &protocols::elscripts::MPI_Slave::available_mem, this );

    master_comm_ = protocols::wum2::EndPointSP( new protocols::wum2::MPI_EndPoint( world, ref_available_mem ) );

		lregister_MPI_Slave(lstate_);
		luabind::globals(lstate_)["slave"] = this;
		luabind::globals(lstate_)["rank"] = world_.rank();
}

void MPI_Slave::go(){
	using namespace utility::lua;

	// create functors to call back functions for automated endpoint processing
  boost::function<void ( protocols::wum2::StatusResponse, int )> ref_listen_wu_sendrecv = boost::bind(  &protocols::wum2::EndPoint::listen_wu_sendrecv, master_comm_, _1, _2);

	// entering main loop
  while( 1 ) {
    master_comm_->check_and_act_clearcommand();
    master_comm_->check_and_act_status_request( ref_listen_wu_sendrecv );

    // check if slave doesn't have enough free memory to run WU
    // (because wu will generate objects that take memory)
    while( mem_limit_ - current_mem() < 2 * reserved_mem_ ) {
      // can't run another WU, otherwise we have to throw away WU from inbuf
      // and we should guarantee all WU sent to slave are run
      
      // only way to free mem is to get rid of outbound stuff
      master_comm_->check_and_act_status_request( ref_listen_wu_sendrecv );
      master_comm_->cleanup_reqs();

      // hopefully we rarely get to here
    }

    // run the first wu in the queue
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
					"tmp_run_env = {} -- delete tmp_env after calling run_on_slave()\n"
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

    // moves inbuf to inq
    // cleans up outbuf
    master_comm_->cleanup_reqs();
  }
}


} //elscripts
} //protocols
#endif
