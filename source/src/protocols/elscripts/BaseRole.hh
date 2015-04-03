// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/BaseRole.hh
/// @brief  BaseRole, which handles a lot of common functions between all the elscripts roles
/// 				including a lot of lua stuff
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_BaseRole_hh
#define INCLUDED_protocols_elscripts_BaseRole_hh

#ifdef USELUA
// this is useless without lua
#include <boost/cstdint.hpp>
#include <boost/unordered_map.hpp>

#include <utility/lua/LuaObject.hh>

#include <core/io/serialization/PipeMap.fwd.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/inputter/InputterStream.fwd.hh>
#include <protocols/inputter/Inputter.fwd.hh>
#include <protocols/outputter/Outputter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/wum2/WorkUnit.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>


namespace protocols {
namespace elscripts {

void lregister_BaseRole( lua_State * lstate );

class BaseRole {
  public:
		BaseRole( boost::uint64_t mem_limit, boost::uint64_t reserved_mem, boost::uint64_t reserved_mem_multiplier):
			mem_limit_(mem_limit),
			reserved_mem_(reserved_mem),
			reserved_mem_multiplier_(reserved_mem_multiplier),
			mover_cache_mem_(0),
			mover_cache_length_(0)
			{
				mover_cache_ = protocols::moves::MoverCacheSP ( new protocols::moves::MoverCache() );
			}
    ~BaseRole(){}
    virtual void go() = 0;

		// memory free after setting some aside for buffer (reserved*multiplier)
    virtual boost::uint64_t available_mem() = 0;

		void reparse_def( std::string const & type, std::string const & name );

  protected:

		// real memory free, ignoring buffer
    virtual boost::uint64_t current_mem()=0;

    boost::uint64_t mover_cache_mem() { return mover_cache_mem_; }

    void update_mover_cache_mem();

		// opens new lua state, exports definition tables, registers classes
		void lua_init();
		void instantiate_inputters();
		void instantiate_inputterstream();
		void instantiate_output();
		void instantiate_movers();
		void instantiate_filters();
		void instantiate_scorefxns();
		void instantiate_tasks();
		void instantiate_workunits();
		void register_calculators();


	protected:
		lua_State * lstate_;
    boost::uint64_t mem_limit_;
    boost::uint64_t reserved_mem_;
    boost::uint64_t reserved_mem_multiplier_;

    protocols::inputter::InputterStreamSP inputterstream_;

    protocols::moves::MoverCacheSP mover_cache_;
    // holds number of members of cache, so that when it changes we can update cache size
    // pathetic logic
    int mover_cache_length_; 
    boost::uint64_t mover_cache_mem_;

    core::io::serialization::PipeMapSP structure_cache_;

		// these are just wrappers around lua tables
		// the lua tables contain boost shared ptrs to cpp allocated mem
		// so when the table entry is set to to nil, lua GC will delete the boost shared ptr
		// which will then delete the cpp allocated mem
		// this also allows direct manipulation of the table entries from lua
		utility::lua::LuaObject movers_;
		utility::lua::LuaObject filters_;
		utility::lua::LuaObject inputters_;
		utility::lua::LuaObject outputters_;
		utility::lua::LuaObject scorefxns_;
		utility::lua::LuaObject tasks_;
};

} //elscripts
} //protocols
#endif
#endif
