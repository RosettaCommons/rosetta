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
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <test/protocols/sewing/extra_functions.hh>
#include <protocols/sewing/movers/AssemblyMover.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

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



