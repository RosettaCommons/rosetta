// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/Outputter.cc
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a file
/// @author Ken Jung

// Unit Headers
#include <protocols/outputter/Outputter.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace outputter {

static thread_local basic::Tracer TR( "protocols.outputter.Outputter" );

#ifdef USELUA
void lregister_Outputter( lua_State * lstate ) {
	using namespace core::io::serialization;
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("outputter")
		[
			luabind::class_<Outputter>("Outputter")
						.def("write", (void (Outputter::*) (PipeMap&)) &Outputter::write)
						.def("write", (void (Outputter::*) (Pipe&)) &Outputter::write)
						.def("write", (void (Outputter::*) (Pose&)) &Outputter::write)
		]
	];
}
#endif
Outputter::Outputter(){}
Outputter::~Outputter(){}

} // outputter
} // protocols
