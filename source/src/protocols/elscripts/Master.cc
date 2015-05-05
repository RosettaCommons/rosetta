// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/Master.cc
/// @brief  The Master role in elscripts, handles trajectories, generating workunits, processing of results
/// @author Ken Jung

#ifdef USELUA
#include <protocols/elscripts/Master.hh>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>

#include <protocols/wum2/WorkUnit.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/els.OptionKeys.gen.hh>

#include <utility/lua/LuaIterator.hh>

#include <utility/Factory.hh>
#include <protocols/inputter/InputterStream.hh>
#include <protocols/inputter/Inputter.hh>
#include <protocols/outputter/Outputter.hh>

#include <basic/Tracer.hh>

// elscripts master
namespace protocols {
namespace elscripts {

static thread_local basic::Tracer TR( "protocols.elscripts.Master" );

void lregister_Master( lua_State * lstate ) {
	lregister_BaseRole( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<Master, BaseRole>("Master")
					.def("make_wu_until_limit", &Master::make_wu_until_limit)
					.def("make_wu", (void (Master::*) ( std::string const &, core::pose::PoseSP, protocols::moves::SerializableStateSP ) ) &Master::make_wu)
					.def("make_wu", (void (Master::*) ( std::string const &, core::io::serialization::PipeSP, protocols::moves::SerializableStateSP ) ) &Master::make_wu)
					.def("make_wu", (void (Master::*) ( std::string const &, core::io::serialization::PipeMapSP, protocols::moves::SerializableStateSP ) ) &Master::make_wu)
					.def("make_wu", (void (*) ( Master *, std::string const &, core::pose::PoseSP ) ) &master_make_wu_nostate)
					.def("make_wu", (void (*) ( Master *, std::string const &, core::io::serialization::PipeSP ) ) &master_make_wu_nostate)
					.def("make_wu", (void (*) ( Master *, std::string const &, core::io::serialization::PipeMapSP ) ) &master_make_wu_nostate)
					.def("interpreter", &Master::interpreter)
					.def("end_traj", &Master::end_traj)
		]
	];
}

Master::Master( int num_trajectories, boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier) :
  num_trajectories_(num_trajectories),
	trajectories_mem_(0),
  num_trajectories_finished_(0),
	mpicounter_(10000000),
	traj_idx_(0),
	last_generate_initial_wu_time_( boost::posix_time::microsec_clock::universal_time() - boost::posix_time::minutes(100)),
	BaseRole( mem_limit, reserved_mem, reserved_mem_multiplier) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// endpoint uses this function to get role-wide memory usage
    boost::function< boost::uint64_t ()> ref_available_mem = boost::bind( &protocols::elscripts::Master::available_mem, this );

    slave_comm_ = protocols::wum2::EndPointSP( new protocols::wum2::EndPoint( ref_available_mem ));


		// initializing inputterstream
		inputterstream_.reset ( new protocols::inputter::InputterStream (
					inputter_rank(),
					// num of masters
					1
					) );

		lua_init();
		lregister_Master(lstate_);
		luabind::globals(lstate_)["master"] = this;
		luabind::globals(lstate_)["num_trajectories"] = num_trajectories_;
		instantiate_tasks();
		instantiate_inputters();
		instantiate_inputterstream();
		instantiate_output();

		register_calculators();
		instantiate_scorefxns();
		instantiate_filters();
		instantiate_movers();
		instantiate_workunits();

		// initializing trajectories and wu counters
    trajectories_ = core::io::serialization::PipeSP( new core::io::serialization::Pipe() );
    for( int i = 0; i < num_trajectories_; i++ ) {
      trajectories_->push_back( core::pose::PoseSP() );
    }

		// setting up a new env for each traj
		std::string action = "DELIM(\n"
			"traj_env = {}\n"
			"for i=0,num_trajectories do\n"
			"	traj_env[i] = {}\n"
			"	setmetatable( traj_env[i], { __index = _G } )\n"
			"	local _ENV = traj_env[i]\n"
			"	traj_idx = i\n"
			"	wu_made = {}\n"
			"	wu_done = {}\n"
			"end\n"
			")DELIM";
    int err = luaL_dostring ( lstate_, action.c_str() );
    if( err == 1) {
      TR << "Setting up trajectory environments failed. Error is:" << std::endl;
      TR << lua_tostring(lstate_, -1) << std::endl;
      std::exit(9);
    }

		luabind::globals(lstate_)["trajectories"] = trajectories_;
}

void Master::interpreter() {
	TR << "Switching to lua interpreter mode!" << std::endl;
	std::string line = "";
	int err = 0;
	while( 1 ) {
		std::cout << "> " << std::flush;
		std::getline(std::cin, line);
		if( line == "quit" ) break;
    err = luaL_dostring ( lstate_, line.c_str() );
    if( err == 1) {
			std::cout << lua_tostring(lstate_, -1) << std::endl;
    }
	}
	TR << "Leaving to lua interpreter mode!" << std::endl;
}

void Master::go(){
	using namespace utility::lua;

	// no while loop since Single Node will be switching between Master::go() and Slave::go()

	mpicounter_++;

	fill_trajectories();

	// dunno if this is needed, int check faster than boost time check?
	if( mpicounter_ >= 10000000 ) {
		// hardcoded, only try and generate initial workunits every 60 seconds
		if ( ( boost::posix_time::microsec_clock::universal_time() - last_generate_initial_wu_time_) > boost::posix_time::seconds( 60 )){

			// generate initial wu
			// user must be aware of memory limits if they use their own function
			int m = luaL_dostring ( lstate_, "loop_every()" );
			if( m == 1) {
				TR << "calling lua function loop_every() failed. Error is:" << std::endl;
				TR << lua_tostring(lstate_, -1) << std::endl;
				std::exit(9);
			}
			last_generate_initial_wu_time_ = boost::posix_time::microsec_clock::universal_time();
		}
		mpicounter_ = 0;
	}

	// process slave results
	if( ! slave_comm_->inq().empty() ) {
		protocols::wum2::WorkUnitSP wu = slave_comm_->inq().pop_front();
		protocols::wum2::WorkUnit_ElScriptsSP castattempt = boost::dynamic_pointer_cast<protocols::wum2::WorkUnit_ElScripts> (wu);
		if( castattempt != 0 ) {
			traj_idx_ = castattempt->trajectory_idx();
			std::string wuname = castattempt->name();

			// export new variables to lua for use in the lua fxn
			luabind::globals(lstate_)["traj_env"][traj_idx_]["pipemap"] = castattempt->pipemap().lock();
			luabind::globals(lstate_)["traj_env"][traj_idx_]["state"] = castattempt->state().lock();

			// calling proceed fxn
			// also increment wu_done here to save overhead of calling into lua vm
			std::string action = "(\n"
				"do\n"
				"	local _ENV = traj_env[)" + boost::lexical_cast<std::string>(traj_idx_) + "(]\n"
				"	if wu_done.)"+wuname+"( == nil then\n"
				"		wu_done.)"+wuname+"( = 0\n"
				"	end\n"
				"	wu_done.)"+wuname+"( = wu_done.)"+wuname+"( + 1\n"
				"	els_setenv(_ENV)\n"
				"	els.workunits.)"+wuname+"DELIM(.proceed_on_master()\n"
				"end\n"
				")DELIM";
			int err = luaL_dostring ( lstate_, action.c_str() );
			if( err == 1) {
				TR << "Calling lua function for workunit " << wuname << " proceed_on_master fxn failed. Error is:" << std::endl;
				TR << lua_tostring(lstate_, -1) << std::endl;
				std::exit(9);
			}
			luabind::globals(lstate_)["traj_env"][traj_idx_]["pipemap"] = luabind::nil;
			luabind::globals(lstate_)["traj_env"][traj_idx_]["state"] = luabind::nil;
			lua_gc(lstate_, LUA_GCCOLLECT, 0);
		}
	}

}

// this uses the "global" traj_idx_
// traj_idx is set before this fxn is called
// this gets rid of traj_idx in lua script
void Master::make_wu( std::string const & wuname, core::pose::PoseSP p, protocols::moves::SerializableStateSP state ) {
	using namespace core::io::serialization;
	PipeSP pipe ( new Pipe() );
	// have to copy the pose here
	// originally, did not deep copy pose for speed reasons
	// but then the problem becomes when boost serialize sends the pose to the slave
	// 2 workunits can point to the same pose!  because thats how it was on the master
	// boost serialize too smart for its own good
	pipe->push_back( core::pose::PoseSP( new core::pose::Pose(*p)) );
	PipeMapSP pmap( new PipeMap() );
	(*pmap)["input"] = pipe;
	if (! state)
		state = protocols::moves::SerializableStateSP(new protocols::moves::SerializableState() );
	make_wu_copied( wuname, pmap, state );
}

void Master::make_wu( std::string const & wuname, core::io::serialization::PipeSP p, protocols::moves::SerializableStateSP state ) {
	using namespace core::io::serialization;
	// have to deep copy, see make_wu( pose ) for reason
	PipeMapSP pmap( new PipeMap() );
	(*pmap)["input"] = core::io::serialization::clone(p);
	if (! state)
		state = protocols::moves::SerializableStateSP(new protocols::moves::SerializableState() );
	make_wu_copied( wuname, pmap, state );
}

void Master::make_wu( std::string const & wuname, core::io::serialization::PipeMapSP pmap, protocols::moves::SerializableStateSP state ) {
	using namespace core::io::serialization;
	// have to deep copy, see make_wu( pose ) for reason
	PipeMapSP newpmap = core::io::serialization::clone(pmap);
	if (! state)
		state = protocols::moves::SerializableStateSP(new protocols::moves::SerializableState() );
	make_wu_copied( wuname, pmap, state );
}

void Master::make_wu_copied( std::string const & wuname, core::io::serialization::PipeMapSP pmap, protocols::moves::SerializableStateSP state) {
	// if you're calling this, pmap should already be a unique instance (deep copied)
	// and state isn't null ^^
	using namespace core::io::serialization;

	protocols::wum2::WorkUnitSP tmp (new protocols::wum2::WorkUnit_ElScripts(
				1, traj_idx_, pmap, state, wuname
				));

	slave_comm_->outq().push_back( tmp );
	if( !luabind::globals(lstate_)["traj_env"][traj_idx_]["wu_made"][wuname] ) {
		luabind::globals(lstate_)["traj_env"][traj_idx_]["wu_made"][wuname] = 0;
	}
	luabind::globals(lstate_)["traj_env"][traj_idx_]["wu_made"][wuname] = luabind::object_cast<int>(luabind::globals(lstate_)["traj_env"][traj_idx_]["wu_made"][wuname]) + 1;
}

void Master::make_wu_until_limit( std::string const & wuname, int num ) {
	using namespace core::io::serialization;
  for( int i = 0; i < num_trajectories_; i++ ){
		if( ! (*trajectories_)[i] ) continue;
		if( !luabind::globals(lstate_)["traj_env"][i]["wu_made"][wuname] ) {
			luabind::globals(lstate_)["traj_env"][i]["wu_made"][wuname] = 0;
		}
		if( luabind::object_cast<int>(luabind::globals(lstate_)["traj_env"][i]["wu_made"][wuname]) < num ) {
			if( mem_limit_ - current_mem() > 2 * reserved_mem_ ) {
				PipeMapSP pmap( new PipeMap() );
				PipeSP pipe ( new Pipe() );
				pipe->push_back( (*trajectories_)[i] );
				(*pmap)["input"] = pipe;
				protocols::moves::SerializableStateSP state( new protocols::moves::SerializableState() );

				protocols::wum2::WorkUnitSP tmp (new protocols::wum2::WorkUnit_ElScripts(
							1, i, pmap, state, wuname
							));
				slave_comm_->outq().push_back( tmp );
				luabind::globals(lstate_)["traj_env"][i]["wu_made"][wuname] = luabind::object_cast<int>(luabind::globals(lstate_)["traj_env"][i]["wu_made"][wuname]) + 1;
			} else {
				return;
			}
		}
	}
}

void Master::end_traj() {
	(*trajectories_)[traj_idx_].reset();
	num_trajectories_finished_++;
	TR << "Finished " << num_trajectories_finished_ << " trajectories." << std::endl;

	// reset the traj env
	std::string action = "(\n"
		"i = )" + boost::lexical_cast<std::string>(traj_idx_) + "DELIM(\n"
		"traj_env[i] = {}\n"
		"setmetatable( traj_env[i], { __index = _G } )\n"
		"do\n"
		"	local _ENV = traj_env[i]\n"
		"	traj_idx = i\n"
		"	wu_made = {}\n"
		"	wu_done = {}\n"
		"end\n"
		")DELIM";
		int err = luaL_dostring ( lstate_, action.c_str() );
	if( err == 1) {
		TR << "Cleaning up trajectory " << traj_idx_ << " environment failed. Error is:" << std::endl;
		TR << lua_tostring(lstate_, -1) << std::endl;
		std::exit(9);
	}
}

void Master::fill_trajectories() {
  std::vector< int > needs_replace;
  for( int i = 0; i < num_trajectories_; i++ ){
    if( ! (*trajectories_)[i] )
      needs_replace.push_back( i );
  }

  while( needs_replace.size() != 0 && inputterstream_->has_pose() ) {
    (*trajectories_)[ needs_replace.back() ] = inputterstream_->get_pose();
    needs_replace.pop_back();
  }
  //request_pool_structures( needs_replace );
  if( needs_replace.size() == num_trajectories_ ) {
		TR << "Elscripts finished successfully" << std::endl;
    exit(9);
	}
}

void Master::update_trajectories_mem(){
#ifdef SERIALIZATION
  std::stringstream s;
//  core::io::serialization::toBinary(s, trajectories_);
  trajectories_mem_ =  s.str().length();
#else
  TR << "Memory usage tracked only if compiled against boost::serialize" << std::endl;
#endif
}

void master_make_wu_nostate( Master * master, std::string const & wuname, core::pose::PoseSP p ) { master->make_wu( wuname, p ); }
void master_make_wu_nostate( Master * master, std::string const & wuname, core::io::serialization::PipeSP p ) { master->make_wu( wuname, p ); }
void master_make_wu_nostate( Master * master, std::string const & wuname, core::io::serialization::PipeMapSP p ){ master->make_wu( wuname, p ); }

} //elscripts
} //protocols
#endif
