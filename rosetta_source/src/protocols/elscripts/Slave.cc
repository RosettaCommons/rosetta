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

static basic::Tracer TR("protocols.elscripts.Slave");

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
				LuaObject dworkunits( luabind::globals(lstate_)["elscripts"]["dworkunits"]);

				luabind::globals(lstate_)["pipemap"] = castattempt->pipemap().lock();
				luabind::globals(lstate_)["state"] = castattempt->state().lock();
				try { 
					luabind::call_function<void>( dworkunits[ castattempt->name() ]["run"].raw() );
				} catch (std::exception & e) {
					TR << "calling lua function for workunit " << castattempt->name() << " run fxn failed failed. Vague error is:" << std::endl;
					TR << e.what() << std::endl;
					TR << lua_tostring(lstate_, -1) << std::endl;
					std::exit(9);
				}
				luabind::globals(lstate_)["pipemap"] = luabind::nil;
				luabind::globals(lstate_)["state"] = luabind::nil;
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
