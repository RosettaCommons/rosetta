// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @brief  test suite for core::chemical::Elements.hh
/// @author Steven Combs


// Test headers
#include <cxxtest/TestSuite.h>
#include <core/chemical/Elements.hh>


using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class ElementsTests : public CxxTest::TestSuite {

public:

	ElementsTests() {}

	virtual ~ElementsTests() {}

	static ElementsTests *createSuite() {
		return new ElementsTests();
	}

	static void destroySuite( ElementsTests *suite ) {
		delete suite;
	}


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. The fixture above
	// gets constructed once before all the tests in this test suite are run.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_elements_to_name() {
		std::string C("C");
		std::string N("N");
		std::string P("P");
		std::string Am("Am");

		TS_ASSERT_EQUALS(core::chemical::element::name_from_elements(element::C), C);
		TS_ASSERT_EQUALS(core::chemical::element::name_from_elements(element::N), N);
		TS_ASSERT_EQUALS(core::chemical::element::name_from_elements(element::P), P);
		TS_ASSERT_EQUALS(core::chemical::element::name_from_elements(element::Am), Am);
	}

	void test_name_to_elements() {
		std::string Pr("Pr");
		std::string Ti("Ti");
		std::string H("H");

		TS_ASSERT_EQUALS( core::chemical::element::Pr, core::chemical::element::elements_from_name(Pr) );
		TS_ASSERT_EQUALS( core::chemical::element::Ti, core::chemical::element::elements_from_name(Ti) );
		TS_ASSERT_EQUALS( core::chemical::element::H, core::chemical::element::elements_from_name(H) );

	}


};

