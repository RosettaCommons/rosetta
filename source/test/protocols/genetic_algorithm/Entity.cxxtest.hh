// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/BumpGrid.cxxtest.hh
/// @brief  test suite for protocols::match::SixDHasher
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/genetic_algorithm/Entity.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// Test util headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <utility/string_util.hh>

#include <sstream>

//Auto Headers
#include <utility/vector1.hh>


using namespace protocols::genetic_algorithm;
using namespace core;


class DummyEntityElementCreator : public EntityElementCreator {
public:
	virtual ~DummyEntityElementCreator() {}
	virtual std::string widget_name() const { return "DUMMY"; }
	virtual EntityElementOP new_entity( std::string const & word );

};


class DummyEntityElement : public EntityElement
{
public:
	typedef EntityElement parent;
public:

	DummyEntityElement() : parent(), dummy_( 0 ) {}
	DummyEntityElement( Size index, Size dummy ) : parent( index ), dummy_( dummy ) {}
	DummyEntityElement( std::string word ) : parent( word )
	{
		std::istringstream iss( word );
		iss >> dummy_;
	}

	virtual ~DummyEntityElement() {}

	virtual EntityElementOP clone() { return EntityElementOP( new DummyEntityElement(*this) );}
	virtual EntityElementOP fresh_instance() { return EntityElementOP( new DummyEntityElement ); }

	virtual Size hash() const { return dummy_; }

	virtual bool operator <  ( EntityElement const & rhs ) const {
		if ( parent::operator < ( rhs ) ) return true;
		if ( parent::operator == ( rhs ) ) return false;

		DummyEntityElement const * rhs_dummy_ptr = dynamic_cast< DummyEntityElement const * > ( &rhs );
		if ( !rhs_dummy_ptr ) {
			utility_exit_with_message( "Dynamic cast to DummyEntityElement failed in operator < with rhs name = " + rhs.name() );
		}
		return dummy_ < rhs_dummy_ptr->dummy_;
	}

	virtual bool operator == ( EntityElement const & rhs ) const
	{
		if ( ! parent::operator == ( rhs ) ) return false;

		DummyEntityElement const * rhs_dummy_ptr = dynamic_cast< DummyEntityElement const * > ( &rhs );
		if ( !rhs_dummy_ptr ) {
			utility_exit_with_message( "Dynamic cast to DummyEntityElement failed in operator == with rhs name = " + rhs.name() );
		}
		return dummy_ == rhs_dummy_ptr->dummy_;
	}

	virtual EntityElement & operator =  ( EntityElement const & rhs ) {
		parent::operator = ( rhs );
		if ( this != &rhs ) {
			DummyEntityElement const * rhs_dummy_ptr = dynamic_cast< DummyEntityElement const * > ( &rhs );
			if ( !rhs_dummy_ptr ) {
				utility_exit_with_message( "Dynamic cast to DummyEntityElement failed in operator = with rhs name = " + rhs.name() );
			}
			dummy_ = rhs_dummy_ptr->dummy_;
		}
		return *this;
	}

	virtual std::string to_string() const {
		return parent::to_string() + utility::to_string( dummy_ );
	}

	virtual std::string name() const {
		DummyEntityElementCreator deec;
		return deec.widget_name();
	}

	Size dummy() const {
		return dummy_;
	}

private:
	Size dummy_;

};


EntityElementOP DummyEntityElementCreator::new_entity( std::string const & word )
{
	return EntityElementOP( new DummyEntityElement( word ) );
}

EntityElementRegistrator< DummyEntityElementCreator > dummyEE_registrator;


// --------------- Test Class --------------- //

class EntityTests : public CxxTest::TestSuite {

public:


	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_DummyElement_ctor()
	{
		DummyEntityElement d;
		TS_ASSERT( d.name() == "DUMMY" );
	}

	Entity
	create_dummy_entity() {
		Entity e;
		utility::vector1< EntityElementOP > traits;
		traits.reserve(10);
		for ( Size ii = 1; ii <= 10; ++ii ) {
			traits.push_back( EntityElementOP( new DummyEntityElement( ii, 2 * ii ) ));
		}
		e.set_traits( traits );
		return e;
	}

	void
	compare_to_gold( Entity const & e ) {
		for ( Size ii = 1; ii <= 10; ++ii ) {
			TS_ASSERT( dynamic_cast< DummyEntityElement const * > ( e.traits()[ ii ].get() ) );
			if ( dynamic_cast< DummyEntityElement const * > ( e.traits()[ ii ].get() ) ) {
				TS_ASSERT( e.traits()[ ii ]->index() == ii );
				TS_ASSERT( (dynamic_cast< DummyEntityElement const & > (*e.traits()[ ii ]) ).dummy() == 2*ii );
			}
		}
	}

	void test_Entity_with_dummy_elements()
	{
		Entity e = create_dummy_entity();
		compare_to_gold( e );
	}

	void test_Entity_output() {
		Entity e = create_dummy_entity();
		std::ostringstream oss;
		oss << e;
		std::string entity_string = oss.str();
		std::string gold_string = "Entity with traits: DUMMY:1:2 DUMMY:2:4 DUMMY:3:6 DUMMY:4:8 DUMMY:5:10 DUMMY:6:12 DUMMY:7:14 DUMMY:8:16 DUMMY:9:18 DUMMY:10:20 and fitness  0.000";
		TS_ASSERT( entity_string == gold_string );
	}

	void test_Entity_from_string() {
		std::string gold_string = "traits DUMMY:1:2 DUMMY:2:4 DUMMY:3:6 DUMMY:4:8 DUMMY:5:10 DUMMY:6:12 DUMMY:7:14 DUMMY:8:16 DUMMY:9:18 DUMMY:10:20 fitness  0.000";
		Entity e( gold_string );
		compare_to_gold( e );
	}

	void test_Entity_equality() {
		Entity e = create_dummy_entity();
		Entity e2 = create_dummy_entity();
		TS_ASSERT( e == e2 );
	}

	void test_Entity_assignment_operator() {
		Entity e = create_dummy_entity();
		Entity e2;
		e2 = e;
		TS_ASSERT( e == e2 );
		e2.traits()[ 1 ]->index( 1234 ); // if the assignment operator performed a deep copy, then e is still gold
		TS_ASSERT( e2.traits()[ 1 ]->index() == 1234 );
		compare_to_gold( e );

	}

	void test_Entity_copy_constructor() {
		Entity e = create_dummy_entity();
		Entity e2(e);
		TS_ASSERT( e == e2 );
		e2.traits()[ 1 ]->index( 1234 ); // if the copy constructor performed a deep copy, then e is still gold
		TS_ASSERT( e2.traits()[ 1 ]->index() == 1234 );
		compare_to_gold( e );
	}

	void test_Entity_fitness_valid_logic() {
		Entity e = create_dummy_entity();
		TS_ASSERT( ! e.fitness_valid() );
		e.set_fitness( 14 );
		TS_ASSERT( e.fitness_valid() );
	}

	void test_Entity_fitness_lt_logic() {
		Entity e = create_dummy_entity();
		Entity e2( e );
		Entity e3( e );
		e.set_fitness( 14 );
		e2.set_fitness( 15 );
		e3.set_fitness( 14 );
		TS_ASSERT( e < e2 );
		TS_ASSERT( ! (e2 < e) );
		TS_ASSERT( ! (e < e3) );
		TS_ASSERT( ! (e3 < e) );
	}

};


