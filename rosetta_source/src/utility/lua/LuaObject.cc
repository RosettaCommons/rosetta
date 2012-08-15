// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=1 noet:
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


#ifdef USELUA
#include <utility/lua/LuaIterator.hh>
#include <utility/lua/LuaObject.hh>


namespace utility {
namespace lua {

LuaObject::operator bool() {
		return object_.is_valid() && luabind::type( object_ ) != LUA_TNIL;
}

LuaIterator LuaObject::begin() const {
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaIterator( luabind::iterator(object_) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		exit(9);
}

LuaObject LuaObject::operator[] ( const std::string & str ) const {
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ str ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		exit(9);
}

// look into luaL_getn
int LuaObject::size() const {
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				int counter = 0;
				for(luabind::iterator i(object_), end; i != end; i++ ){
						counter++;
				}
				return counter;
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		exit(9);
}

LuaObject LuaObject::operator[] ( const char * str ) const {
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ str ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		exit(9);
}

LuaObject LuaObject::operator[] ( int i ) const {
		if (i == 0 ) {
				std::cerr << "----------ERROR---------------"  << std::endl
				<< "Lua tables start indexing from 1, 0 will always return nil" << std::endl;
				exit(9);
		}
		if( luabind::type( object_ ) == LUA_TTABLE ) {
				return LuaObject( luabind::object(object_[ i ]) );
		}
		std::cerr << "----------ERROR---------------"  << std::endl
		<< "Attempted to access index of non-table lua object" << std::endl;
		exit(9);
}



} //lua
} //utility
#endif
