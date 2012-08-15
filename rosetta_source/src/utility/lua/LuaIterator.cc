// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=1 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/lua/LuaIterator.cc
/// @brief A wrapper around a luabind::iterator
/// has skey() and ikey() functions that return a string or int conversion of the key
// 
// this is really just a convenience class
/// @author Ken Jung

#ifdef USELUA
#include <utility/lua/LuaIterator.hh>
#include <utility/lua/LuaObject.hh>

#include <boost/mpl/assert.hpp>

namespace utility {
namespace lua {

// skey and ikey go through LuaObject because it has more verbose error handling of failed cast
std::string LuaIterator::skey() {
		return LuaObject( iterator_.key() ).to<std::string>();
}

int LuaIterator::ikey() {
		return LuaObject( iterator_.key() ).to<int>();
}

luabind::iterator LuaIterator::raw() {
		return iterator_;
}

bool LuaIterator::operator==(LuaIterator & other) {
		return iterator_ == other.raw();
}

bool LuaIterator::operator!=(LuaIterator & other) {
		return iterator_ != other.raw();
}

LuaObject *  LuaIterator::operator -> () {
		// wish i could somehow detect -> usage at compile time
		std::cerr << "-------- ERROR ---------" << std::endl
				<< "\t -> not supported by LuaIterator, use (*itr). instead" << std::endl;
		exit(9);
		// will never get past here
		LuaObject * tmp = new LuaObject();
		return tmp;
}
LuaObject LuaIterator::operator * () {
		return LuaObject(*iterator_);
}

LuaIterator LuaIterator::operator++(int) {
		LuaIterator tmp( iterator_ );
		iterator_++;
		return tmp;
}

LuaIterator & LuaIterator::operator++() {
		iterator_++;
		return *this;
}



} //lua
} //utility
#endif
