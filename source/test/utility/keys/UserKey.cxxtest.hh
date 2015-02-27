// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/keys/UserKey.cxxtest.hh
/// @brief  UserKey unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/UserKey.hh>

#include <platform/types.hh>


// --- set up some helper classes, in their own namespace
namespace userkey {

class X
{};


class XKey : public utility::keys::UserKey< X >
{

	friend class X;

	typedef  utility::keys::UserKey< X >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	explicit
	XKey(
		std::string const & id_a = std::string(),
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{}


	/// @brief Index constructor
	inline
	XKey(
		Size const & index_a,
		std::string const & id_a = std::string(),
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( index_a, id_a, identifier_a, code_a )
	{}


	/// @brief Clone this
	inline
	XKey *
	clone() const
	{
		return new XKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~XKey()
	{}

};

} // namespace userkey


// grant access to our helper classes
using namespace userkey;

class UserKeyTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_UserKey_general() {

		XKey::Size n_key( 0 );

		// make first key
		XKey key1( ++n_key, "Key1" );
		TS_ASSERT( key1.id() == "Key1" );

		// make second key
		XKey key2( ++n_key, "Key2" );
		TS_ASSERT( key2.id() == "Key2" );

		// compare keys
		TS_ASSERT( key1 < key2 );
	}

};


