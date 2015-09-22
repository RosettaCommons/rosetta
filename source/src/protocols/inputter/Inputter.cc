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

// Unit Headers
#include <protocols/inputter/Inputter.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace inputter {

static THREAD_LOCAL basic::Tracer TR( "protocols.inputter.Inputter" );

#ifdef USELUA
void lregister_Inputter( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("inputter")
		[
			luabind::class_<Inputter>("Inputter")
		]
	];
}
#endif
Inputter::Inputter(){
	offset_ = false;
}
Inputter::~Inputter(){}

} // inputter
} // protocols
