// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/sewing/movers/SewingAssemblyMoverTests.cxxtest.hh
/// @brief  Tests for the sewing AssemblyMover class
/// @author Minnie Langlois (minnie@email.unc.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/sewing/movers/AssemblyMover.hh>


// Core Headers

// Protocol Headers
#include <basic/Tracer.hh>

#include <core/init_util.hh> // AUTO IWYU For core_init

using namespace protocols::sewing;

static basic::Tracer TR("SewingAssemblyMoverTests");


class SewingAssemblyMoverTests : public CxxTest::TestSuite {
	//Define Variables
private:
	movers::AssemblyMover a_mover_;

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}







};



