// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/InputterStream.cc
/// @brief InputterStream holds a list of of streams and does VERY basic things with them
/// like controls whether structures are duplicated across masters
/// and controls whether to take the list sequentially or round robin
/// @author Ken Jung

// Unit Headers
#include <protocols/inputter/InputterStream.hh>
#include <protocols/inputter/Inputter.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace inputter {

static thread_local basic::Tracer TR( "protocols.inputter.InputterStream" );

#ifdef USELUA
void lregister_InputterStream( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("inputter")
		[
			luabind::class_<InputterStream>("InputterStream")
		]
	];
}
#endif

InputterStream::~InputterStream(){}

bool InputterStream::has_pose() {
	input_itr itr = inputters_.begin();
	while ( itr != inputters_.end() ) {
		if ( (**itr).offset() ) {
			if ( (**itr).has_nth_pose( num_masters_ ) ) {
				return true;
			}
		} else if ( (**itr).has_nth_pose( master_rank_ ) ) {
			return true;
		}
		// if we got here, this inputter is all out of poses for us
		itr = inputters_.erase(itr);
	}
	return false;
}

core::pose::PoseSP InputterStream::get_pose() {
	// it is assumed that you call has_pose before this
	// so the first inputter always is valid
	if ( inputters_.front()->offset() ) {
		return inputters_.front()->get_nth_pose( num_masters_ );
	} else {
		return inputters_.front()->get_nth_pose( master_rank_ );
	}
	// should never reach here
	return core::pose::PoseSP();
}

void InputterStream::add_inputter( InputterSP inputter ) {
	inputters_.push_back( inputter );
}

void InputterStream::parse_def( utility::lua::LuaObject const & /*def*/ ) {
}

} // inputter
} // protocols
