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

#ifdef USEBOOSTSERIALIZE
#include <boost/archive/binary_oarchive.hpp>
#endif

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

static basic::Tracer TR("protocols.elscripts.Master");

void lregister_Master( lua_State * lstate ) {
	lregister_BaseRole( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<Master, BaseRole>("Master")
					.def("make_wu_until_limit", &Master::make_wu_until_limit)
					.def("make_wu", &Master::make_wu)
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
		instantiate_tasks();
		instantiate_inputters();
		instantiate_inputterstream();
		instantiate_output();

		register_calculators();
		instantiate_scorefxns();
		instantiate_filters();
		instantiate_movers();

		// initializing trajectories and wu counters
    trajectories_ = core::io::serialization::PipeSP( new core::io::serialization::Pipe() );
		luabind::globals(lstate_)["elscripts"]["wumade"] = luabind::newtable( lstate_ );
		luabind::globals(lstate_)["elscripts"]["wumade"] = luabind::newtable( lstate_ );
		luabind::globals(lstate_)["elscripts"]["wufinished"] = luabind::newtable( lstate_ );
		wumade_.raw( luabind::globals(lstate_)["elscripts"]["wumade"] );
		wufinished_.raw( luabind::globals(lstate_)["elscripts"]["wufinished"] );
    for( int i = 0; i < num_trajectories_; i++ ) {
      trajectories_->push_back( core::pose::PoseSP() );
			luabind::globals(lstate_)["elscripts"]["wumade"][i] = luabind::newtable( lstate_ );
			luabind::globals(lstate_)["elscripts"]["wufinished"][i] = luabind::newtable( lstate_ );
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
			int m = luaL_dostring ( lstate_, "generate_initial_wus()" );
			if( m == 1) {
				TR << "calling lua function generate_initial_wus() failed. Error is:" << std::endl;
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
			if( ! wufinished_[traj_idx_][wuname] ) {
				luabind::globals(lstate_)["elscripts"]["wufinished"][traj_idx_][wuname] = 0;
			}
			luabind::globals(lstate_)["elscripts"]["wufinished"][traj_idx_][wuname] = wufinished_[traj_idx_][wuname].to<int>() + 1;
			// export new variables to lua for use in the lua fxn
			luabind::globals(lstate_)["traj_idx"] = traj_idx_;
			luabind::globals(lstate_)["pipemap"] = castattempt->pipemap().lock();
				luabind::globals(lstate_)["state"] = castattempt->state().lock();
			luabind::globals(lstate_)["wufinished"] = wufinished_[traj_idx_].raw();
			luabind::globals(lstate_)["wumade"] = wumade_[traj_idx_].raw();

			LuaObject dworkunits( luabind::globals(lstate_)["elscripts"]["dworkunits"]);
			try { 
				luabind::call_function<void>( dworkunits[ wuname]["proceed"].raw() );
			} catch (std::exception & e) {
				TR << "calling lua function for workunit " << wuname << " proceed fxn failed failed. Vague error is:" << std::endl;
				TR << e.what() << std::endl;
				TR << lua_tostring(lstate_, -1) << std::endl;
				lua_pop(lstate_, 1);
				std::exit(9);
			}
			luabind::globals(lstate_)["traj_idx"] = luabind::nil;
			luabind::globals(lstate_)["wufinished"] = luabind::nil;
			luabind::globals(lstate_)["wumade"] = luabind::nil;
			luabind::globals(lstate_)["pipemap"] = luabind::nil;
				luabind::globals(lstate_)["state"] = luabind::nil;
			lua_gc(lstate_, LUA_GCCOLLECT, 0);
		}
	}

}

// this uses the "global" traj_idx_
// traj_idx is set before this fxn is called
// this gets rid of traj_idx in lua script
void Master::make_wu( std::string const & wuname, core::pose::PoseSP p) {
	using namespace core::io::serialization;
	PipeMapSP pmap( new PipeMap() );
	PipeSP pipe ( new Pipe() ); 
	pipe->push_back( core::pose::PoseSP(p) );
	(*pmap)["input"] = pipe;
	protocols::moves::SerializableStateSP state( new protocols::moves::SerializableState() );

	protocols::wum2::WorkUnitSP tmp (new protocols::wum2::WorkUnit_ElScripts(
				1, traj_idx_, pmap, state, wuname
				));

	slave_comm_->outq().push_back( tmp );
	if( ! wumade_[traj_idx_][wuname] ) {
		luabind::globals(lstate_)["elscripts"]["wumade"][traj_idx_][wuname] = 0;
	}
	luabind::globals(lstate_)["elscripts"]["wumade"][traj_idx_][wuname] = wumade_[traj_idx_][wuname].to<int>() + 1;
}

void Master::make_wu_until_limit( std::string const & wuname, int num ) {
	using namespace core::io::serialization;
  for( int i = 0; i < num_trajectories_; i++ ){
		if( ! (*trajectories_)[i] ) continue; 
		if( ! wumade_[i][wuname] ) {
			luabind::globals(lstate_)["elscripts"]["wumade"][i][wuname] = 0;
		}
		while( wumade_[i][ wuname ].to<int>() < num ) {
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
				luabind::globals(lstate_)["elscripts"]["wumade"][i][wuname] = wumade_[i][wuname].to<int>() + 1;
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

	luabind::globals(lstate_)["elscripts"]["wumade"][traj_idx_] = luabind::newtable( lstate_ );
	luabind::globals(lstate_)["elscripts"]["wufinished"][traj_idx_] = luabind::newtable( lstate_ );
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
#ifdef USEBOOSTSERIALIZE
  std::stringstream s;
  boost::archive::binary_oarchive oa(s);
  oa << trajectories_;
  trajectories_mem_ =  s.str().length();
#else
  TR << "Memory usage tracked only if compiled against boost::serialize" << std::endl;
#endif
}


} //elscripts
} //protocols
#endif
