// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/components/DivideAndConquerorTests.cxxtest.hh
/// @brief  Test suite for breaking up a denovo structure and folding the pieces
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/components/DivideAndConqueror.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

static basic::Tracer TR("DivideAndConquerorTests");

class DivideAndConquerorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}


	void test_rossman() {
		using protocols::denovo_design::components::BuildPhases;
		using protocols::denovo_design::components::DivideAndConqueror;
		using protocols::denovo_design::components::StructureData;
		using protocols::denovo_design::components::StructureDataFactory;

		std::stringstream xml;
		xml << "<StructureData name=\"rossman\" length=\"97\" pose_length=\"99\" >" << std::endl
			<< "<ResidueRange name=\"E3\" start=\"2\" group=1 ss=\"LEEEEEE\" abego=\"XBBBBBB\" nterm=\"0\" cterm=\"1\" upper_segment=\"E3_H1\" />" << std::endl
			<< "<ResidueRange name=\"E3_H1\" start=\"8\" group=2 ss=\"LLL\" abego=\"XXX\" nterm=\"1\" cterm=\"1\" lower_segment=\"E3\" upper_segment=\"H1\" />" << std::endl
			<< "<ResidueRange name=\"H1\" start=\"11\" group=3 ss=\"HHHHHHHHHHHHHH\" abego=\"AAAAAAAAAAAAAA\" nterm=\"1\" cterm=\"1\" lower_segment=\"E3_H1\" upper_segment=\"H1_E4\" />" << std::endl
			<< "<ResidueRange name=\"H1_E4\" start=\"25\" group=4 ss=\"LL\" abego=\"XX\" nterm=\"1\" cterm=\"1\" lower_segment=\"H1\" upper_segment=\"E4\" />" << std::endl
			<< "<ResidueRange name=\"E4\" start=\"27\" group=1 ss=\"EEEEE\" abego=\"BBBBB\" nterm=\"1\" cterm=\"1\" lower_segment=\"H1_E4\" upper_segment=\"E4_H2\" />" << std::endl
			<< "<ResidueRange name=\"E4_H2\" start=\"32\" group=5 ss=\"LL\" abego=\"XX\" nterm=\"1\" cterm=\"1\" lower_segment=\"E4\" upper_segment=\"H2\" />" << std::endl
			<< "<ResidueRange name=\"H2\" start=\"34\" group=6 ss=\"HHHHHHHHHHHHHH\" abego=\"AAAAAAAAAAAAAA\" nterm=\"1\" cterm=\"1\" lower_segment=\"E4_H2\" upper_segment=\"H2_E2\" />" << std::endl
			<< "<ResidueRange name=\"H2_E2\" start=\"48\" group=7 ss=\"LLLL\" abego=\"XXXX\" nterm=\"1\" cterm=\"1\" lower_segment=\"H2\" upper_segment=\"E2\" />" << std::endl
			<< "<ResidueRange name=\"E2\" start=\"52\" group=1 ss=\"EEEEEE\" abego=\"BBBBBB\" nterm=\"1\" cterm=\"1\" lower_segment=\"H2_E2\" upper_segment=\"E2_H3\" />" << std::endl
			<< "<ResidueRange name=\"E2_H3\" start=\"58\" group=8 ss=\"LL\" abego=\"XX\" nterm=\"1\" cterm=\"1\" lower_segment=\"E2\" upper_segment=\"H3\" />" << std::endl
			<< "<ResidueRange name=\"H3\" start=\"60\" group=9 ss=\"HHHHHHHHHHHHHH\" abego=\"AAAAAAAAAAAAAA\" nterm=\"1\" cterm=\"1\" lower_segment=\"E2_H3\" upper_segment=\"H3_E1\" />" << std::endl
			<< "<ResidueRange name=\"H3_E1\" start=\"74\" group=10 ss=\"LL\" abego=\"XX\" nterm=\"1\" cterm=\"1\" lower_segment=\"H3\" upper_segment=\"E1\" />" << std::endl
			<< "<ResidueRange name=\"E1\" start=\"76\" group=1 ss=\"EEEEEE\" abego=\"BBBBBB\" nterm=\"1\" cterm=\"1\" lower_segment=\"H3_E1\" upper_segment=\"E1_H4\" />" << std::endl
			<< "<ResidueRange name=\"E1_H4\" start=\"82\" group=11 ss=\"LLL\" abego=\"XXX\" nterm=\"1\" cterm=\"1\" lower_segment=\"E1\" upper_segment=\"H4\" />" << std::endl
			<< "<ResidueRange name=\"H4\" start=\"85\" group=12 ss=\"HHHHHHHHHHHHHHL\" abego=\"AAAAAAAAAAAAAAX\" nterm=\"1\" cterm=\"0\" lower_segment=\"E1_H4\" />" << std::endl
			<< "<Pairing type=\"StrandPairing\" orient1=\"1\" orient2=\"1\" shift=\"0\" segments=\"E1,E2\" />" << std::endl
			<< "<Pairing type=\"StrandPairing\" orient1=\"1\" orient2=\"1\" shift=\"0\" segments=\"E2,E3\" />" << std::endl
			<< "<Pairing type=\"StrandPairing\" orient1=\"1\" orient2=\"1\" shift=\"1\" segments=\"E3,E4\" />" << std::endl
			<< "</StructureData>";

		StructureData const sd = StructureDataFactory::get_instance()->create_from_xml( xml );
		TR << sd << std::endl;
		DivideAndConqueror mr_conq;
		BuildPhases const st_segments = mr_conq.divide_and_conquer( sd );

		// 4-phase build expected
		TS_ASSERT_EQUALS( st_segments.size(), 4 );

		// Phase 1 : strands 3-4, helix 1
		BuildPhases::const_iterator seg = st_segments.begin();
		TS_ASSERT_EQUALS( *seg, boost::assign::list_of ("E3")("E3_H1")("H1")("H1_E4")("E4") );

		// Phase 2 : add helix 2, strand 2
		++seg;
		TS_ASSERT_EQUALS( *seg, boost::assign::list_of ("E4")("E4_H2")("H2")("H2_E2")("E2") );

		// Phase 3 : add helix 3, strand 1
		++seg;
		TS_ASSERT_EQUALS( *seg, boost::assign::list_of ("E2")("E2_H3")("H3")("H3_E1")("E1") );

		++seg;
		TS_ASSERT_EQUALS( *seg, boost::assign::list_of ("E1")("E1_H4")("H4") );
	}

	void test_ntf2()
	{
		using protocols::denovo_design::components::BuildPhases;
		using protocols::denovo_design::components::DivideAndConqueror;
		using protocols::denovo_design::components::StructureData;
		using protocols::denovo_design::components::StructureDataFactory;

		std::ifstream xml_in( "protocols/denovo_design/components/DivideAndConquerorTests_NTF2.xml" );
		StructureData const sd = StructureDataFactory::get_instance()->create_from_xml( xml_in );
		TR << sd << std::endl;

		DivideAndConqueror conq;
		conq.set_start_segments( boost::assign::list_of("E04") );
		BuildPhases const st_segments = conq.divide_and_conquer( sd );
		TR << sd << std::endl;
		TR << st_segments << std::endl;
	}


};



