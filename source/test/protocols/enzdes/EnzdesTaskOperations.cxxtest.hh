// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/enzdestaskoperations.cxxtest.hh
/// @brief  test suite for enzdestaskoperations
/// @details
/// @author Steve Bertolani sjbertolani@ucdavis.edu

// Test headers
#include <cxxtest/TestSuite.h>

// initialization headers
#include <test/core/init_util.hh>

//basic
#include <basic/Tracer.hh>

// EnzDestaskOperations
#include <protocols/enzdes/EnzdesTaskOperations.hh>

static basic::Tracer TR("test.protocols.enzdes.EnzdesTaskOperations");

class EnzdesTaskOperations_Tests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_ProteinLigandInterfaceUpweighter_gettersetter() {
		/*
		yes, this is a tad bit silly....  but then again, before this test the get_weight returned a bool....
		*/
		//Creating ProteinLigandInterfaceUpweighter
		protocols::enzdes::ProteinLigandInterfaceUpweighter * protliginterf_ = new protocols::enzdes::ProteinLigandInterfaceUpweighter;

		// check initialization to 1.0
		TS_ASSERT( protliginterf_->get_weight() == 1.0 );

		// now check set weight works
		protliginterf_->set_weight( 2.0 );
		TS_ASSERT( protliginterf_->get_weight() == 2.0);
	}
}; // class
