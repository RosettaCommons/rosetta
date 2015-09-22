// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/MPI_Master.cc
/// @brief  The MPI_Master role in elscripts, handles trajectories, generating workunits, processing of results
/// @author Ken Jung

#if defined (USEBOOSTMPI) && defined (USELUA)
#include <protocols/elscripts/MPI_Master.hh>

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

static THREAD_LOCAL basic::Tracer TR( "protocols.elscripts.MPI_Master" );

void lregister_MPI_Master( lua_State * lstate ) {
	lregister_Master( lstate );
	luabind::module(lstate, "protocols")
		[
		luabind::namespace_("elscripts")
		[
		luabind::class_<MPI_Master, Master>("MPI_Master")
		]
		];
}

MPI_Master::MPI_Master( boost::mpi::communicator world, std::vector<int> slaves, int num_trajectories, boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier) :
	world_(world),
	last_status_sweep_time_( boost::posix_time::microsec_clock::universal_time() - boost::posix_time::minutes(100)),
	slaves_( slaves ),
	Master( num_trajectories, mem_limit, reserved_mem, reserved_mem_multiplier) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	// endpoint uses this function to get role-wide memory usage
	boost::function< boost::uint64_t ()> ref_available_mem = boost::bind( &protocols::elscripts::MPI_Master::available_mem, this );

	slave_comm_ = protocols::wum2::EndPointSP( new protocols::wum2::MPI_EndPoint( world, ref_available_mem ));

	// initializing inputterstream
	inputterstream_.reset ( new protocols::inputter::InputterStream (
		inputter_rank(),
		// num of masters
		(option[OptionKeys::els::num_traj]() / option[OptionKeys::els::traj_per_master]() ) +
		(!( option[OptionKeys::els::num_traj]() % option[OptionKeys::els::traj_per_master]() == 0 ))
		) );

	lregister_MPI_Master(lstate_);
	luabind::globals(lstate_)["master"] = this;
	luabind::globals(lstate_)["rank"] = world_.rank();

	instantiate_inputterstream();
}

void MPI_Master::go(){
	using namespace utility::lua;

	// create functors to call back functions for automated endpoint processing
	boost::function<void ( protocols::wum2::StatusResponse, int )> ref_slave_listen_wu_sendrecv = boost::bind(  &protocols::wum2::EndPoint::listen_wu_sendrecv, slave_comm_, _1, _2);
	boost::function<bool ( protocols::wum2::StatusResponse )> ref_slave_initiate_wu_sendrecv = boost::bind(  &protocols::wum2::EndPoint::initiate_wu_sendrecv, slave_comm_, _1);
	//  boost::function<void ( protocols::wum2::StatusResponse, int )> ref_pool_listen_wu_sendrecv = boost::bind(  &protocols::wum2::EndPoint::listen_wu_sendrecv, pool_comm_, _1, _2);
	// boost::function<bool ( protocols::wum2::StatusResponse )> ref_pool_initiate_wu_sendrecv = boost::bind(  &protocols::wum2::EndPoint::initiate_wu_sendrecv, pool_comm_, _1);

	// setting intial sweep time such that the first sweep will always occur immediately

	// entering main loop
	while ( 1 ) {
		mpicounter_++;

		// pool_comm_->check_and_act_clearcommand();
		//pool_comm_->check_and_act_status_request( f );

		// only send status requests as often as the shortest_wu would take
		if ( ( boost::posix_time::microsec_clock::universal_time() - last_status_sweep_time_) > boost::posix_time::seconds(basic::options::option[basic::options::OptionKeys::els::shortest_wu]()) ) {
			for ( int i = 0; i < slaves_.size(); i++ ) {
				if ( ! slave_comm_->has_open_status( slaves_[i] ) ) {
					// only send status requests to slaves that have responded
					// if they havent responded yet, spamming them with more messages won't help
					// probably should retry a few times after a timeout
					// removing slave from slave list is unecessary, but we do lose workunits sent there
					slave_comm_->send_status_request( slaves_[i] );
				}
			}
			last_status_sweep_time_ = boost::posix_time::microsec_clock::universal_time();
		}

		fill_trajectories();

		if ( mpicounter_ >= 10000000 ) {
			// hardcoded, only try and generate initial workunits every 60 seconds
			if ( ( boost::posix_time::microsec_clock::universal_time() - last_generate_initial_wu_time_) > boost::posix_time::seconds( 60 ) ) {

				// generate initial wu
				// user must be aware of memory limits if they use their own function
				int m = luaL_dostring ( lstate_, "loop_every()" );
				if ( m == 1 ) {
					TR << "calling lua function loop_every() failed. Error is:" << std::endl;
					TR << lua_tostring(lstate_, -1) << std::endl;
					std::exit(9);
				}
				last_generate_initial_wu_time_ = boost::posix_time::microsec_clock::universal_time();
			}
			mpicounter_ = 0;
		}

		//  send/recv WU from slaves
		slave_comm_->act_on_status_response( ref_slave_initiate_wu_sendrecv );

		//  send/recv WU from pool
		//pool_comm_->check_and_act_status_request( ref_pool_listen_wu_sendrecv );

		// process slave results
		protocols::wum2::WorkUnitSP wu = slave_comm_->inq().pop_front();
		protocols::wum2::WorkUnit_ElScriptsSP castattempt = boost::dynamic_pointer_cast<protocols::wum2::WorkUnit_ElScripts> (wu);
		if ( castattempt != 0 ) {
			traj_idx_ = castattempt->trajectory_idx();
			std::string wuname = castattempt->name();
			// export new variables to lua for use in the lua fxn
			luabind::globals(lstate_)["traj_env"][traj_idx_]["pipemap"] = castattempt->pipemap().lock();
			luabind::globals(lstate_)["traj_env"][traj_idx_]["state"] = castattempt->state().lock();

			// calling proceed fxn
			// also increment wu_done here to save overhead of calling into lua vm
			std::string action = "(\n"
				"do\n"
				"\tlocal _ENV = traj_env[)" + boost::lexical_cast<std::string>(traj_idx_) + "(]\n"
				"\tif wu_done.)"+wuname+"( == nil then\n"
				"\t\twu_done.)"+wuname+"( = 0\n"
				"\tend\n"
				"\twu_done.)"+wuname+"( = wu_done.)"+wuname+"( + 1\n"
				"\tels_setenv(_ENV)\n"
				"\tels.workunits.)"+wuname+"DELIM(.proceed_on_master()\n"
				"end\n"
				")DELIM";
			int err = luaL_dostring ( lstate_, action.c_str() );
			if ( err == 1 ) {
				TR << "Calling lua function for workunit " << wuname << " proceed_on_master fxn failed. Error is:" << std::endl;
				TR << lua_tostring(lstate_, -1) << std::endl;
				std::exit(9);
			}
			luabind::globals(lstate_)["traj_env"][traj_idx_]["pipemap"] = luabind::nil;
			luabind::globals(lstate_)["traj_env"][traj_idx_]["state"] = luabind::nil;
			lua_gc(lstate_, LUA_GCCOLLECT, 0);
		}

		// process pool results
		/*WorkUnitSP wu = pool_comm_->inq_popfront();
		WorkUnit_RequestStructSP castattempt = boost::dynamic_pointer_cast<WorkUnit_RequestStruct> (wu);
		if( castattempt != 0 ) {
		trajectories[cast_attempt->traj_idx()] = cast_attempt->replacement_pose();
		}
		*/

		// cleans up any completed mpi::reqs and their buffers
		slave_comm_->cleanup_reqs();
		//pool_comm_->cleanup_reqs();
	}
}

/*void MPI_Master::request_pool_structures( std::vector < int > needs_replace ) {
for( int i = 0; i < needs_replace.size(); i++ ) {
request_pool_structure( i );
}
}

void MPI_Master::request_pool_structure( int traj_idx ) {
protocols::wum2::WorkUnitSP request_struct_wu( new protocols::wum2::WorkUnit_RequestStruct() );
pool_comm_.push_back( request_struct_wu );
}*/

} //elscripts
} //protocols
#endif
