// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/SilentFileInputter.hh
/// @brief An inputter that takes a list of silent files
/// @author Ken Jung

#ifndef INCLUDED_protocols_inputter_SilentFileInputter_hh
#define INCLUDED_protocols_inputter_SilentFileInputter_hh

// Unit Headers
#include <protocols/inputter/SilentFileInputter.fwd.hh>
#include <protocols/inputter/Inputter.hh>
#include <core/io/silent/SilentFileData.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <deque>


namespace protocols {
namespace inputter {

#ifdef USELUA
void lregister_SilentFileInputter( lua_State * lstate );
#endif

class SilentFileInputter : public Inputter {

public:
	SilentFileInputter():multiplier_(1), multiply_over_all_(true), curr_idx_(0){}
	virtual ~SilentFileInputter();

	// throw away n-1 poses and return the nth one
	// necessary to prevent duplication of input across different masters
	// of course, default is 1 for non-mpi scenarios
	core::pose::PoseSP get_nth_pose( int n=1 );
	bool has_nth_pose( int n=1 );

#ifdef USELUA
		// need to pass in a map of the previous inputters, as inputters will call inputters
		// not a huge fan of this, but implementing it correctly and more extensibly will be a lot more work
		// and fpd doesn't think there will be much customizing requested in inputters
		void parse_def( utility::lua::LuaObject const & def,
						utility::lua::LuaObject const & tasks,
						utility::lua::LuaObject & inputters );
		virtual void lregister( lua_State * lstate );
#endif

	// factory functions
	InputterSP create();
	static std::string name() {
		return "SilentFileInputter";
	}


private:
	std::deque< std::pair< int, std::string> > tags_;
	int multiplier_;
	bool multiply_over_all_;
	core::Size curr_idx_;
	core::io::silent::SilentFileData sfd_;

}; // end SilentFileInputter

} // inputter
} // protocols


#endif //INCLUDED_protocols_inputter_Inputter_HH
