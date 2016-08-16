// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/AutoKey.cxxtest.hh
/// @brief  AutoKey unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/AutoKey.hh>


// --- set up some helper classes, in their own namespace
namespace autokey {

class X_
{};


class XKey_ : public ::utility::keys::AutoKey< X_ >
{

	friend class X_;
	typedef  ::utility::keys::AutoKey< X_ >  Super;

public: // Creation


	/// @brief Default constructor
	inline
	explicit
	XKey_(
		std::string const & id_a = std::string(),
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{}


	/// @brief Clone this
	inline
	XKey_ *
	clone() const
	{
		return new XKey_( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~XKey_()
	{}


};

} // namespace autokey


// grant access to our helper classes
using namespace autokey;

class AutoKeyTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_AutoKey_general() {

		// make first key
		XKey_ key1( "Key1" );
		TS_ASSERT( key1.id() == "Key1" );

		// make second key
		XKey_ key2( "Key2" );
		TS_ASSERT( key2.id() == "Key2" );

		// compare keys
		TS_ASSERT( key1 < key2 );
	}

};


