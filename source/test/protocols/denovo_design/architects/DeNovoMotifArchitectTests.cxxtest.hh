// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/components/MotifArchitectBase.cxxtest.hh
/// @brief  Tests for StructureArchitect and MotifArchitectBase classes
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

// Core Headers
#include <core/pose/Pose.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <iterator>

static THREAD_LOCAL basic::Tracer TR("MotifArchitectBase");

using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::components;

class DeNovoMotifArchitectTests : public CxxTest::TestSuite {
public:
	void setUp()
	{
		core_init();
	}

	void tearDown()
	{}


	void test_apply()
	{
	}

	void test_parse_my_tag()
	{
	}

	void test_motif_base()
	{
		std::string const test_id = "test_id3";

		DeNovoMotifArchitect dummy( test_id );
		TS_ASSERT_EQUALS( std::distance( dummy.motifs_begin(), dummy.motifs_end() ), 0 );

		//This will induce a compile error -- I still need to write this and I don't wnat to forget
	}

};

