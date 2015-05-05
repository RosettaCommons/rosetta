// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/serialization/PipeMap.fwd.hh
/// @brief Quick typedef of a map of strings to PipeSP
/// @author Ken Jung

#ifndef INCLUDED_core_io_serialization_PipeMap_fwd_hh
#define INCLUDED_core_io_serialization_PipeMap_fwd_hh

#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#endif

#include <core/pose/Pose.fwd.hh>
#include <core/io/serialization/Pipe.fwd.hh>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <map>

namespace core {
namespace io {
namespace serialization {

typedef std::map< std::string, PipeSP > PipeMap;
typedef boost::shared_ptr< PipeMap > PipeMapSP;
typedef boost::weak_ptr< PipeMap > PipeMapWP;

#ifdef USELUA
void lregister_PipeMap( lua_State * lstate );
#endif
void insert( PipeMapSP p, std::string const & pipename, core::pose::PoseSP pose );
PipeSP at( PipeMapSP p, std::string const & pipename );
PipeMapSP clone( PipeMapSP p);
PipeMapSP inputPipeMap( core::pose::PoseSP p );

} //serialization
} //io
} //core

#endif


