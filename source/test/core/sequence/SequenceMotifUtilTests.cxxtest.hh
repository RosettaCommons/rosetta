// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/sequence/SequenceMotifUtilTests.cxxtest.hh
/// @brief  Test utility functions for sequence motifs.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/sequence/sequence_motif.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SequenceMotifUtilTests");


class SequenceMotifUtilTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_split_sequence_motif(){

		std::string motif = "N[^P][ST]VTR[%POLAR]";

		utility::vector1< std::string > motifSP = core::sequence::split_sequence_motif( motif );

		TS_ASSERT_EQUALS( motifSP[1], "N+" );
		TS_ASSERT_EQUALS( motifSP[2], "^P");
		TS_ASSERT_EQUALS( motifSP[3], "ST");
		TS_ASSERT_EQUALS( motifSP[4], "+VTR+");
		TS_ASSERT_EQUALS( motifSP[5], "%POLAR");

	}






};
