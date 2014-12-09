// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/lua/LuaObject.cc
/// @brief A wrapper around a luabind::object 
/// has bool conversion, [] support, nested table support, and a .to<T>() conversion function
/// since a luabind::object is actually a pointer to the stack of a luastate object,
/// if the luastate object is destroyed, before this wrapper class is destroyed, i don't know what will happen
// 
// this is really just a convenience class
/// @author Ken Jung


#include <utility/lua/LuaIterator.hh>
#include <utility/lua/LuaObject.hh>
#include <utility/exit.hh>


namespace utility {
namespace lua {

LuaObject::operator bool() {
#ifdef USELUA
		return object_.is_valid() && luabind::type( object_ ) != LUA_TNIL;
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return false;
#endif
}

LuaIterator LuaObject::begin() const {
#ifdef USELUA
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaIterator( luabind::iterator(object_) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		utility_exit_with_message("");
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return LuaIterator();
#endif
}

LuaObject LuaObject::operator[] ( const std::string & /*str*/ ) const {
#ifdef USELUA
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ str ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		utility_exit_with_message("");
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return LuaObject();
#endif
}

// look into luaL_getn
int LuaObject::size() const {
#ifdef USELUA
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				int counter = 0;
				for(luabind::iterator i(object_), end; i != end; i++ ){
						counter++;
				}
				return counter;
		}
		std::cerr << "----------ERROR---------------"  << std::endl;
		utility_exit_with_message("Attempted to access index of non-table lua object" );
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return 0;
#endif
}

LuaObject LuaObject::operator[] ( const char * /*str*/ ) const {
#ifdef USELUA
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ str ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl;
		utility_exit_with_message("Attempted to access index of non-table lua object" );
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return LuaObject();
#endif
}

LuaObject LuaObject::operator[] ( int /*i*/ ) const {
#ifdef USELUA
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ i ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		utility_exit_with_message("");
#else
		utility_exit_with_message("Can't use LuaObject without compiling with USELUA flag" );
		return LuaObject();
#endif
}



} //lua
} //utility
