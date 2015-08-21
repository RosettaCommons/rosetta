// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/InputterStream.hh
/// @brief InputterStream holds a list of of streams and does VERY basic things with them
/// like controls whether structures are duplicated across masters
/// and controls whether to take the list sequentially or round robin
/// @author Ken Jung

#ifndef INCLUDED_protocols_inputter_InputterStream_hh
#define INCLUDED_protocols_inputter_InputterStream_hh

// Unit Headers
#include <protocols/inputter/InputterStream.fwd.hh>
#include <protocols/inputter/Inputter.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/lua/LuaObject.hh>
#include <list>

namespace protocols {
namespace inputter {

#ifdef USELUA
void lregister_InputterStream( lua_State * lstate );
#endif

class InputterStream {
	typedef std::list< InputterSP >::iterator input_itr;

public:
	InputterStream( int master_rank, int num_masters) :
		master_rank_(master_rank),
		num_masters_(num_masters){}

	virtual ~InputterStream();

	virtual bool has_pose();
	virtual core::pose::PoseSP get_pose();
	virtual void add_inputter( InputterSP inputter );
	virtual void parse_def( utility::lua::LuaObject const & def );
	virtual int size() { return inputters_.size(); }

private:
	std::list< InputterSP > inputters_;
	int master_rank_;
	int num_masters_;

}; // end class

} // inputter
} // protocols


#endif //INCLUDED_protocols_inputter_Inputter_HH
