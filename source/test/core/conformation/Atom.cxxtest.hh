// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreTest.cxxtest.hh
/// @brief  unified scoring test.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>


#include <core/types.hh>

// Unit headers
#include <core/conformation/Atom.hh>

#include <test/UTracer.hh>

//Auto Headers


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.ConformationAtom.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace conformation;

///////////////////////////////////////////////////////////////////////////
/// @name EnergyMapTest
/// @brief: Test the functionality of the EnergyMap class
///////////////////////////////////////////////////////////////////////////
class ConformationAtomTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	/// @brief test that the default constructor sets the atom type to 0
	/// and positions the atom at the origin.
	void test_ConformationAtom_default_constructor() {
		Atom atom;
		TS_ASSERT_EQUALS( atom.type(), 0 );
		TS_ASSERT_EQUALS( atom.mm_type(), 0 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().y(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().z(), 0.0 );
	}

	/// @brief test that the atomtype constructor initializes
	/// the type to the input value and the coordinates to the origin
	void test_ConformationAtom_atomtype_constructor() {
		Atom atom( 1, 1 );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().y(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().z(), 0.0 );
	}

	/// @brief test that the copy constructor copies both the type
	/// and the coordinate data
	void test_ConformationAtom_copy_constructor() {
		Atom atom( 1, 1 );
		Vector point( 0.5, -0.5, 2.25 );
		atom.xyz( point );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.5 );
		TS_ASSERT_EQUALS( atom.xyz().y(), -0.5 );
		TS_ASSERT_EQUALS( atom.xyz().z(), 2.25 );

		Atom copy_atom( atom );
		TS_ASSERT_EQUALS( copy_atom.type(), 1 );
		TS_ASSERT_EQUALS( copy_atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( copy_atom.xyz().x(), 0.5 );
		TS_ASSERT_EQUALS( copy_atom.xyz().y(), -0.5 );
		TS_ASSERT_EQUALS( copy_atom.xyz().z(), 2.25 );
	}

	/// @brief test that atom returns the correct type
	void test_ConformationAtom_type_getter() {
		Atom atom( 1, 1 );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
	}

	/// @brief test that the atom type setter method stores the correct type
	void test_ConformationAtom_type_setter() {
		Atom atom( 1, 1 );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		atom.type( 2 );
		TS_ASSERT_EQUALS( atom.type(), 2 );
		atom.mm_type( 2 );
		TS_ASSERT_EQUALS( atom.mm_type(), 2 );
	}

	/// @brief test that the xyz getter method returns the correct coordinates.
	/// If this test fails than almost all other tests will fail.
	void test_ConformationAtom_xyz_getter() {
		Atom atom( 1, 1 );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().y(), 0.0);
		TS_ASSERT_EQUALS( atom.xyz().z(), 0.0 );
	}

	/// @brief test that the xyz setter method correctly repositions the atom.
	/// If this test fails, then the copy constructor test will also fail.
	void test_ConformationAtom_xyz_setter() {
		Atom atom( 1, 1 );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.0 );
		TS_ASSERT_EQUALS( atom.xyz().y(), 0.0);
		TS_ASSERT_EQUALS( atom.xyz().z(), 0.0 );

		Vector point( 0.5, -0.5, 2.25 );
		atom.xyz( point );
		TS_ASSERT_EQUALS( atom.type(), 1 );
		TS_ASSERT_EQUALS( atom.mm_type(), 1 );
		TS_ASSERT_EQUALS( atom.xyz().x(), 0.5 );
		TS_ASSERT_EQUALS( atom.xyz().y(), -0.5 );
		TS_ASSERT_EQUALS( atom.xyz().z(), 2.25 );

	}

};
