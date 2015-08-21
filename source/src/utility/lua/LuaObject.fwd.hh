// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/lua/LuaObject.fwd.hh
/// @brief A wrapper around a luabind::object
/// has bool conversion, [] support, nested table support, and a .to<T>() conversion function
/// since a luabind::object is actually a pointer to the stack of a luastate object,
/// if the luastate object is destroyed, before this wrapper class is destroyed, i don't know what will happen
//
// this is really just a convenience class
/// @author Ken Jung

#ifndef INCLUDED_utility_lua_LuaObject_fwd_hh
#define INCLUDED_utility_lua_LuaObject_fwd_hh

namespace utility {
namespace lua {
class LuaObject;
}
}

#endif
