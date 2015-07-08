// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/serialization/PipeMap.cc
/// @brief Quick typedef of a map of strings to PipeSP
/// @author Ken Jung

#include <core/io/serialization/PipeMap.fwd.hh>
#include <core/pose/Pose.hh>

namespace core {
namespace io {
namespace serialization {
#ifdef USELUA
void lregister_PipeMap( lua_State * lstate ) {
	luabind::module(lstate, "core")
	[
		luabind::namespace_("io")
		[
			luabind::namespace_("serialization")
			[
				luabind::class_<PipeMap, PipeMapSP>("PipeMap")
					.def(luabind::constructor<>())
					.def("at", (PipeSP (*) ( PipeMapSP, std::string const & )) &at)
					.def("insert", &insert)
					.def("clone", (PipeMapSP (*) ( PipeMapSP)) &clone),
				luabind::def("inputPipeMap", &inputPipeMap )
			]
		]
	];
}

void insert( PipeMapSP p, std::string const & pipename, core::pose::PoseSP pose ) {
	if (p->find( pipename ) == p->end() )
		(*p)[ pipename ] = PipeSP( new Pipe );
	(*p)[pipename]->push_back( pose );
}

PipeSP at( PipeMapSP p, std::string const & pipename ) {
	return (*p)[pipename];
}
// this is a deep copy, every PoseSP is dereferenced and copied
PipeMapSP clone( PipeMapSP p ) {
	PipeMapSP newpipemap = PipeMapSP( new PipeMap);
	for( PipeMap::iterator itr = p->begin(), end=p->end(); itr != end; ++itr ) {
		(*newpipemap)[itr->first] = PipeSP( new Pipe );
		for( Pipe::iterator jtr = (*itr->second).begin(), jend=(*itr->second).end(); jtr != jend; ++jtr ){
			(*(*newpipemap)[itr->first]).push_back( core::pose::PoseSP( new core::pose::Pose( **jtr ) ) );
		}
	}
	return newpipemap;
}

PipeMapSP inputPipeMap( core::pose::PoseSP p ) {
	PipeMapSP pmap = PipeMapSP( new PipeMap);
	(*pmap)["input"] = PipeSP( new Pipe );
	(*pmap)["input"]->push_back(p);
	return pmap;
}
#endif
} //serialization
} //io
} //core

