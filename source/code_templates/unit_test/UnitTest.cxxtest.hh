// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  --path--/--class--.cxxtest.hh
/// @brief  --brief--
/// @author --name-- (--email--)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("--class--");


class --class-- : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


--test_functions--



};



