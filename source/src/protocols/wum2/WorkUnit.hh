// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum2/WorkUnitBase.hh
/// @brief Base work unit (abstract) for wum2 and some commonly used derived work units
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_WorkUnit_hh
#define INCLUDED_protocols_wum2_WorkUnit_hh

#include <protocols/wum2/WorkUnit.fwd.hh>
#include <core/io/serialization/PipeMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <core/pose/Pose.hh>
#include <core/types.hh>

#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#endif

namespace protocols {
namespace wum2 {

#ifdef USELUA
void lregister_WorkUnit( lua_State * lstate );
void lregister_WorkUnit_Wait( lua_State * lstate );
void lregister_WorkUnit_ElScripts( lua_State * lstate );
#endif

/// @brief  The base class for all work units, this is abstract

class WorkUnit {

public:
	// boost serialize requires empty constructor
	WorkUnit(){}
	// mpi rank of master node, index of trajectory on that master node
	WorkUnit( core::Size master,
			core::Size trajectory_idx );

	virtual ~WorkUnit (){}

	/// @brief Run the workunit
	virtual void run() = 0;

	/// @brief Print WU details to the stream, single line by default
	virtual void print( std::ostream & out, bool verbose = false ) const ;

	/// @brief Set the unixtime of the start of the execution of this WorkUnit
	void set_run_start();

	/// @brief Set the unixtime of the stop of the execution of this WorkUnit
	void set_run_stop();

	/// @brief Returns the difference between unix start and stop times
	core::Size get_run_time();

	// Accessors/Mutators
	void id( int id ) { id_ = id; }
	int id() { return id_; }

	void master( int master ) { master_ = master; }
	int master() { return master_; }

	void trajectory_idx( int trajectory_idx ) { trajectory_idx_ = trajectory_idx; }
	int trajectory_idx() { return trajectory_idx_; }

	void prioritize( bool prioritize ) { prioritize_ = prioritize; }
	bool prioritize() { return prioritize_; }

	// links a cache with the WU so it can use it in run()
	// cache is not seralized nor transferred with WU
	void link_cache( protocols::moves::MoverCacheSP cache) { cache_ = cache; }

protected:

		// whether the WU gets pushed onto the front vs back of the queue
		bool prioritize_;

		// Generates unique id for this WU, from rank and unix timestamp
		void create_unique_id();

		// Unique id for this WU, from rank and unix timestamp
		int id_;

		// For MPI, what master and what trajectory on that master created this WU
		// master = mpi rank of node where WU was created
		int master_;
		int trajectory_idx_;

		/// Important unixtimes
		core::Size  unixtime_creation_;
		core::Size  unixtime_start_;
		core::Size  unixtime_stop_;

		// Cache pointer, is only valid after link_cache() is called
		// is not serialized, is not sent with WU
		protocols::moves::MoverCacheSP cache_;
};  // class WorkUnit


/// @brief WorkUnit that sleeps for X seconds
class WorkUnit_Wait: public WorkUnit {
public:
	// boost serialize requires empty constructor
	WorkUnit_Wait(){}
    WorkUnit_Wait( core::Size master,
			core::Size trajectory_idx,
			long wait_time);

	~WorkUnit_Wait(){}

    void run();

	void wait_time( long wait_time ) { wait_time_ = wait_time; }
	long wait_time() { return wait_time_; }

private:

	// in seconds
	long wait_time_;
};  // class WorkUnit_Wait

// Basic work unit that holds state for ElScripts
class WorkUnit_ElScripts : public WorkUnit {
public:
	// boost serialize requires empty constructor
	WorkUnit_ElScripts(){}
	WorkUnit_ElScripts( core::Size master,
			core::Size trajectory_idx,
			core::io::serialization::PipeMapSP p,
			protocols::moves::SerializableStateSP state,
			std::string name );

	~WorkUnit_ElScripts(){}
    virtual void run();

	// no copy constructor needed, shallow copy is fine for current use

    core::io::serialization::PipeMapWP pipemap() {
		return core::io::serialization::PipeMapWP(pipemap_);
	}

    protocols::moves::SerializableStateWP state() {
		return protocols::moves::SerializableStateWP(state_);
	}

	// cpp doesnt let you overload by return type boohoo
    protocols::moves::SerializableStateCWP const_state() {
		return protocols::moves::SerializableStateCWP(state_);
	}

	// just renames pipes from old name to new name
	// Undefined, commenting out to fix PyRosetta build  void rename_pipes( core::io::serialization::PipeMapSP p, std::map< std::string, std::string > new_names );

	std::string name() { return name_; }
	void name( std::string name ) { name_ = name; }

protected:

	// makes sense to put name here instead of at base class, as this WU is used as the result WU
	// it should hold the name of the WU that generated the result it is holding
	std::string name_;
	core::io::serialization::PipeMapSP pipemap_;
	protocols::moves::SerializableStateSP state_;
};  // class WorkUnit_ElScripts
}
}

#endif
