// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/AtomTests.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Gordon Lemmon


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Element.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Platform Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <ostream>

//Auto Headers


using std::endl;
using std::string;
using basic::Tracer;
using core::Vector;
using core::chemical::AtomICoor;
using core::chemical::Atom;
using core::chemical::ResidueType;

static Tracer TR("core.chemical.AtomTests.cxxtest");

class AtomTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_atom_initialization() {

		Vector xyz(1,2,3);
		core::chemical::ElementOP element;
		Atom atom ("test", "mm_test", 2, element, 1, xyz);

		TS_ASSERT_EQUALS("test", atom.name());
		TS_ASSERT_EQUALS("mm_test", atom.mm_name());
		TS_ASSERT_EQUALS(2, (int) atom.mm_atom_type_index());
		TS_ASSERT_EQUALS(1, atom.charge());
		TS_ASSERT_EQUALS(xyz, atom.ideal_xyz());

	}

	void test_atom_setters(){
		Atom atom;
		atom.name("test");
		atom.mm_name("mm_test");
		atom.atom_type_index(1);
		atom.mm_atom_type_index(2);
		atom.charge(1);
		atom.charge(1);
		atom.ideal_xyz(Vector(1,2,3));

		TS_ASSERT_EQUALS("test", atom.name());
		TS_ASSERT_EQUALS("mm_test", atom.mm_name());
		TS_ASSERT_EQUALS(1, (int) atom.atom_type_index());
		TS_ASSERT_EQUALS(2, (int) atom.mm_atom_type_index());
		TS_ASSERT_EQUALS(1, atom.charge());
		TS_ASSERT_EQUALS(Vector(1,2,3), atom.ideal_xyz());
	}

};
