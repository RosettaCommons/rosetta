// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/Inputter.hh
/// @brief An inputter, a class that returns a poseSP
/// @author Ken Jung

#ifndef INCLUDED_protocols_inputter_Inputter_hh
#define INCLUDED_protocols_inputter_Inputter_hh

// Unit Headers
#include <protocols/inputter/Inputter.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/lua/LuaObject.hh>

#include <string>

namespace protocols {
namespace inputter {

#ifdef USELUA
void lregister_Inputter( lua_State * lstate );
#endif

class Inputter {

public:
	Inputter();
	virtual ~Inputter();

	// throw away n-1 poses and return the nth one
	// necessary to prevent duplication of input across different masters
	// of course, default is 1 for non-mpi scenarios
	virtual core::pose::PoseSP get_nth_pose( int n=1 ) = 0;
	virtual bool has_nth_pose( int n=1 ) = 0;


	// this says whether any poses have been returned from here or not
	// ie whether or not the masterrank offset has been set
	// this is required because if inputter is split across masters
	// 1) get_nth_pose( my_master_rank )
	// 2) get_nth_pose( num_masters) for every pose after that
	// so need to know when 1) has been done
	bool offset() {
		return offset_;
	}

#ifdef USELUA
		// need to pass in a map of the previous inputters, as inputters will call inputters
		// not a huge fan of this, but implementing it correctly and more extensibly will be a lot more work
		// and fpd doesn't think there will be much customizing requested in inputters
		virtual void parse_def( utility::lua::LuaObject const & def,
						utility::lua::LuaObject const & tasks,
						utility::lua::LuaObject & inputters ) = 0;	
		virtual void lregister( lua_State * lstate )=0;
#endif

	// factory functions
	virtual InputterSP create() = 0;
	static std::string name() {
		return "UNDEFINED NAME";
	}


protected:
	bool offset_;

}; // end Inputter base class

} // inputter
} // protocols


#endif //INCLUDED_protocols_inputter_Inputter_HH
