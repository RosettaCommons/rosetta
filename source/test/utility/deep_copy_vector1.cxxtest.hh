// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/deep_copy_vector1.cxxtest.hh
/// @brief  Unit tests for utility::deep_copy_vector1
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Headers
#include <iostream> // Put this first to allow debug print-outs in project headers
#include <set>

// Project headers
#include <cxxtest/TestSuite.h>
#include <utility/deep_copy_vector1.hh>

#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>

// --- set up some helper classes/functions for these tests
// --- (these go in a special namespace to be used from this file only!)

namespace deep_copy_vector1_test {

class ClonableBase;
typedef utility::pointer::shared_ptr< ClonableBase > ClonableBaseOP;
typedef utility::pointer::shared_ptr< ClonableBase const > ClonableBaseCOP;

class ClonableBase {
public:

	ClonableBase( int ident = 0 ):
		ident_( ident )
	{}

	virtual
	~ClonableBase() = default;

	virtual
	ClonableBaseOP clone() const {
		return ClonableBaseOP( new ClonableBase( *this ) );
	}

public: // Data members

	int ident_ = 0;

};

ClonableBaseOP deep_copy( ClonableBase const & source ) { return source.clone(); }


class HolderClass {
public:
	HolderClass() = default;
	HolderClass( core::Size size ) :
		items(size)
	{}

	utility::deep_copy_vector1< ClonableBaseOP > items;
};


} // namespace deep_copy_vector1_test


class DeepCopyVector1Tests : public CxxTest::TestSuite {


public:

	void test_construct() {
		// Just checking that the other vector1 constructors are properly transcluded.
		using namespace deep_copy_vector1_test;
		HolderClass deflt;
		HolderClass size( 2 );
		TS_ASSERT_EQUALS( deflt.items.size(), 0 );
		TS_ASSERT_EQUALS( size.items.size(), 2 );
		TS_ASSERT_EQUALS( size.items[1], nullptr );
		TS_ASSERT_EQUALS( size.items[2], nullptr );
	}

	void test_copying() {
		using namespace deep_copy_vector1_test;

		HolderClass source;
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(1) );
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(2) );
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(3) );

		HolderClass copy( source );

		TS_ASSERT_EQUALS( source.items.size(), 3 );
		TS_ASSERT_EQUALS( copy.items.size(), 3 );

		for ( core::Size ii(1); ii <= 3; ++ii ) {
			TS_ASSERT_EQUALS( copy.items[ii]->ident_, ii );
			TS_ASSERT_DIFFERS( source.items[ii].get(), copy.items[ii].get() ); // Different pointers -> different objects.
		}
	}


	void test_assignment() {
		using namespace deep_copy_vector1_test;

		HolderClass source;
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(1) );
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(2) );
		source.items.push_back( utility::pointer::make_shared< ClonableBase >(3) );

		HolderClass copy;
		copy.items.push_back( utility::pointer::make_shared< ClonableBase >(4) );

		copy = source;

		TS_ASSERT_EQUALS( source.items.size(), 3 );
		TS_ASSERT_EQUALS( copy.items.size(), 3 );

		for ( core::Size ii(1); ii <= 3; ++ii ) {
			TS_ASSERT_EQUALS( copy.items[ii]->ident_, ii );
			TS_ASSERT_DIFFERS( source.items[ii].get(), copy.items[ii].get() ); // Different pointers -> different objects.
		}
	}

	// Return a new HolderClass by value to test move semantics.
	deep_copy_vector1_test::HolderClass make_new_holder( utility::vector1< deep_copy_vector1_test::ClonableBase* > & orig_address) {
		using namespace deep_copy_vector1_test;
		HolderClass retval;
		for ( core::Size ii(1); ii <= 4; ++ii ) {
			ClonableBaseOP entry = utility::pointer::make_shared< ClonableBase >(ii);
			orig_address.push_back( entry.get() );
			retval.items.push_back( entry );
		}
		return retval;
	}

	void test_move_semantics() {
		using namespace deep_copy_vector1_test;
		utility::vector1< ClonableBase* > orig_address1;
		HolderClass move_construct( make_new_holder(orig_address1) );

		TS_ASSERT_EQUALS( move_construct.items.size(), 4 );
		for (  core::Size ii(1); ii <= 4; ++ii ) {
			TS_ASSERT_EQUALS( move_construct.items[ii].get(), orig_address1[ii] );
		}

		utility::vector1< ClonableBase* > orig_address2;
		HolderClass move_assign;
		move_assign = make_new_holder(orig_address2);

		TS_ASSERT_EQUALS( move_assign.items.size(), 4 );
		for (  core::Size ii(1); ii <= 4; ++ii ) {
			TS_ASSERT_EQUALS( move_assign.items[ii].get(), orig_address2[ii] );
		}
	}

};


