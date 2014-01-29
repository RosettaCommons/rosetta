// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/OrbitalTests.cxxtest.hh
/// @brief unit tests for ResidueType
/// @author Gordon Lemmon


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Orbital.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>

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
using core::chemical::ResidueType;
using core::chemical::Orbital;
using core::chemical::Orbital;
using core::chemical::orbitals::ICoorOrbitalData;

static Tracer TR("core.chemical.OrbitalTests.cxxtest");

class OrbitalTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_orbital_initialization() {

		Vector xyz(1,2,3);
		Orbital orbital("test", 1, xyz);

 		TS_ASSERT_EQUALS("test", orbital.name());
 		TS_ASSERT_EQUALS(1, (int) orbital.orbital_type_index());
 		TS_ASSERT_EQUALS(xyz, orbital.ideal_xyz());
 	//	TS_ASSERT_EQUALS(0, orbital.icoor().phi() + orbital.icoor().theta() + orbital.icoor().distance() ); // icoor doesn't have an operator==
 		//TS_ASSERT_EQUALS(0, orbital.new_icoor().phi() + orbital.new_icoor().theta() + orbital.new_icoor().distance() ); // icoor doesn't have an operator==

	}

	void test_orbital_setters(){
		Orbital orbital;
		orbital.name("test");
		orbital.orbital_type_index(1);
		orbital.ideal_xyz(Vector(1,2,3));
		//orbital.icoor(ICoorOrbitalData(1,2,3,'a','b','c'));
	//	orbital.new_icoor(ICoorOrbitalData(4,5,6,'d','e','f'));

 		TS_ASSERT_EQUALS("test", orbital.name());
 		TS_ASSERT_EQUALS(1, (int) orbital.orbital_type_index());
 		///TS_ASSERT_EQUALS(Vector(1,2,3), orbital.ideal_xyz());
 		//TS_ASSERT_EQUALS(6, orbital.icoor().phi() + orbital.icoor().theta() + orbital.icoor().distance() ); // icoor doesn't have an operator==
 		//TS_ASSERT_EQUALS(15, orbital.new_icoor().phi() + orbital.new_icoor().theta() + orbital.new_icoor().distance() ); // icoor doesn't have an operator==
	}

};
