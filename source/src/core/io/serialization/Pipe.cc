// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/serialization/Pipe.cc
/// @brief Quick typedef of a vector of PoseSP
/// @author Ken Jung

#include <core/io/serialization/Pipe.fwd.hh>
#include <core/pose/Pose.hh>

#ifdef USELUA
#include <luabind/iterator_policy.hpp>
#endif

namespace core {
namespace io {
namespace serialization {
// this is a crutch function, luabind needs a wrapper to a stl to return a stl itr
Pipe & each( PipeSP p ) {
	return *p;
}

core::pose::PoseSP at( PipeSP p, int idx ){
	return (*p)[idx];
}

// this is a deep copy, every PoseSP is dereferenced and copied
PipeSP clone( PipeSP p ) {
	PipeSP newpipe = PipeSP( new Pipe);
	for( Pipe::iterator jtr = p->begin(), end = p->end(); jtr != end; ++jtr ){
		newpipe->push_back( core::pose::PoseSP( new core::pose::Pose( **jtr ) ) );
	}
	return newpipe;
}
#ifdef USELUA
void lregister_Pipe( lua_State * lstate ) {
	luabind::module(lstate, "core")
	[
		luabind::namespace_("io")
		[
			luabind::namespace_("serialization")
			[
				luabind::class_<Pipe>("Pipe")
					.def("each", (Pipe & (*) (PipeSP) ) &each, luabind::return_stl_iterator)
					.def("at", (core::pose::PoseSP (*) (PipeSP, int idx) ) &at )
					.def("size", &Pipe::size)
			]
		]
	];
}
#endif
} //serialization
} //io
} //core

