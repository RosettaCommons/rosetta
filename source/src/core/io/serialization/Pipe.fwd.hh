// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/serialization/Pipe.fwd.hh
/// @brief Quick typedef of a vector of PoseSP
/// @author Ken Jung

#ifndef INCLUDED_core_io_serialization_Pipe_fwd_hh
#define INCLUDED_core_io_serialization_Pipe_fwd_hh

#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#endif

#include <boost/shared_ptr.hpp>
#include <core/pose/Pose.fwd.hh>
#include <vector>

namespace core {
namespace io {
namespace serialization {

typedef std::vector< core::pose::PoseSP > Pipe;
typedef boost::shared_ptr< Pipe > PipeSP;

#ifdef USELUA
void lregister_Pipe( lua_State * lstate );
#endif

int size( PipeSP p );
core::pose::PoseSP at( PipeSP p, int idx );
PipeSP clone( PipeSP p );

} //serialization
} //io
} //core

#endif


