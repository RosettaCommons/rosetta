// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/lua/LuaObject.hh
/// @brief A wrapper around a luabind::object
/// has bool conversion, [] support, nested table support, and a .to<T>() conversion function
/// since a luabind::object is actually a pointer to the stack of a luastate object,
/// if the luastate object is destroyed, before this wrapper class is destroyed, i don't know what will happen
//
// this is really just a convenience class
/// @author Ken Jung

#ifndef INCLUDED_utility_lua_LuaObject_hh
#define INCLUDED_utility_lua_LuaObject_hh

#include <utility/lua/LuaIterator.fwd.hh>
#include <utility/lua/LuaObject.fwd.hh>
#include <utility/exit.hh>
#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#endif

#include <iostream>


namespace utility {
namespace lua {

class LuaObject {

public:
#ifdef USELUA
    LuaObject( luabind::object object) :
						object_(object) {}

				luabind::object raw() const {
						return object_;
				}

				void raw(luabind::object object) {
						object_ = object;
				}

#endif

	LuaObject(){}

	~LuaObject(){}

	operator bool();

	LuaIterator begin() const;

	int size() const;

	LuaObject operator[] ( std::string const & str ) const;

	LuaObject operator[] ( const char *  str ) const;

	LuaObject operator[] ( int i ) const;

	template <class T> T to() const {
#ifdef USELUA
						if( ! object_.is_valid() ) {
								std::cerr << "------Error in casting LuaObject to C++ type!-------" << std::endl
										<< "\tThis LuaObject does not refer to a real Lua object." << std::endl
										<< "\tLuaObject must be constructed with an initialized luabind::object." << std::endl
										<< "\tDefault empty constructor for luabind::object is NOT initialized." << std::endl;
								utility_exit_with_message("");
						}
						try {
								return luabind::object_cast<T>(object_);
						}
						catch (std::exception) {
								// first figure out the lua type
								// these are defined in lua.h as LUA_TBOOLEAN, etc
								int luatype = luabind::type( object_ );
								std::string luatypename;
								switch(luatype){
										case 0:
												luatypename = "nil"; break;
										case 1:
												luatypename = "boolean"; break;
										case 2:
												luatypename = "light userdata"; break;
										case 3:
												luatypename = "number"; break;
										case 4:
												luatypename = "string"; break;
										case 5:
												luatypename = "table"; break;
										case 6:
												luatypename = "function"; break;
										case 7:
												luatypename = "userdata"; break;
								}
								std::string cpptypename;
								// hardcode some common types
								if( typeid(T) == typeid(int) ) {
										cpptypename = "int";
								} else if( typeid(T) == typeid(std::string) ) {
										cpptypename = "std::string";
								} else {
										// Yes I know this will result in different output on diff compileres
										// and is usually a mangled name
										// but any debug info is better than none
										cpptypename = typeid(T).name();
								}
								std::cerr << "------Error in casting LuaObject to C++ type!-------" << std::endl
										<< "\tCast from Lua type '" << luatypename
										<< "' to C++ type '" << cpptypename <<  "' failed!" << std::endl;
								utility_exit_with_message("");
						}
						utility_exit_with_message("");
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		// AMW: cppcheck will flag this extraneous return following a utility_exit
		// The trouble is that acknowledging this will kill the windows.cl.PyRosetta build
		T t;
		return t;
#endif
	}

private:
#ifdef USELUA
				luabind::object object_;
#endif

};

} //lua
} //utility
#endif
