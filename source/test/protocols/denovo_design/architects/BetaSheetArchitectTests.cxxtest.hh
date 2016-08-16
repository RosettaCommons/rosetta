// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/BetaSheetArchitectTests.cxxtest.hh
/// @brief  Test suite for BetaSheet Architect
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Protocol Headers
#include <protocols/denovo_design/architects/BetaSheetArchitect.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/SheetDB.hh>

// Core Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR("BetaSheetArchitectTests");

using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::components;

class BetaSheetArchitectTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_simple()
	{
		std::string const sheet_id = "sheet";
		std::string const s1_id = "sheet.s1";
		std::string const s2_id = "sheet.s2";
		std::string const s3_id = "sheet.s3";

		std::stringstream xml;
		xml << "<BetaSheetArchitect name=" << sheet_id << ">" << std::endl;
		xml << "<StrandArchitect name=s1 length=4 orientation=U register_shift=0 />" << std::endl;
		xml << "<StrandArchitect name=s2 length=5 orientation=U register_shift=-1 />" << std::endl;
		xml << "<StrandArchitect name=s3 length=4 orientation=U register_shift=0 />" << std::endl;
		xml << "</BetaSheetArchitect>" << std::endl;
		utility::tag::TagCOP tag = utility::tag::Tag::create( xml );
		basic::datacache::DataMap data;

		BetaSheetArchitect sheet( "" );
		sheet.parse_my_tag( tag, data );

		core::Real random = 0.25001;
		StructureDataOP sd_ptr = sheet.design( core::pose::Pose(), random );
		TS_ASSERT( sd_ptr );
		StructureData const & sd = *sd_ptr;

		TS_ASSERT( sd.has_segment( s1_id ) );
		TS_ASSERT( sd.has_segment( s2_id ) );
		TS_ASSERT( sd.has_segment( s3_id ) );
		TS_ASSERT_EQUALS( sd.segment( s1_id ).ss(), "LEEEEL" );
		TS_ASSERT_EQUALS( sd.segment( s1_id ).abego(), "XBBBBX" );
		TS_ASSERT_EQUALS( sd.segment( s2_id ).ss(), "LEEEEEL" );
		TS_ASSERT_EQUALS( sd.segment( s2_id ).abego(), "XBBBBBX" );
		TS_ASSERT_EQUALS( sd.segment( s3_id ).ss(), "LEEEEL" );
		TS_ASSERT_EQUALS( sd.segment( s3_id ).abego(), "XBBBBX" );
		TS_ASSERT_EQUALS( sd.ss(), "LEEEELLEEEEELLEEEEL" );
		TS_ASSERT_EQUALS( sd.abego(), "XBBBBXXBBBBBXXBBBBX" );
	}

	/* TL: commented until if/when my SheetDB can be fully available in master void nest_with_db()
	{
	using protocols::denovo_design::components::SheetList;

	std::string const sheet_id = "sheet";
	std::string const s1_id = "sheet.s1";
	std::string const s2_id = "sheet.s2";
	std::string const s3_id = "sheet.s3";

	std::stringstream xml;
	xml << "<BetaSheetArchitect name=" << sheet_id << " sheet_db=\"/work/tlinsky/sheet_db/clustered\" >" << std::endl;
	xml << "<Strand name=s1 length=4 orientation=U register_shift=0 />" << std::endl;
	xml << "<Strand name=s2 length=5 orientation=U register_shift=-1 />" << std::endl;
	xml << "<Strand name=s3 length=4 orientation=U register_shift=0 />" << std::endl;
	xml << "</BetaSheetArchitect>" << std::endl;
	utility::tag::TagCOP tag = utility::tag::Tag::create( xml );
	basic::datacache::DataMap data;

	BetaSheetArchitect sheet( "" );
	sheet.parse_my_tag( tag, data );

	// pull this sheet from sheetDB
	protocols::denovo_design::components::SheetDB db;
	db.set_db_path( "/work/tlinsky/sheet_db/clustered" );
	Lengths const test_lengths = boost::assign::list_of (4)(5)(4);
	StrandOrientations const test_orientations = boost::assign::list_of (UP)(UP)(UP);
	RegisterShifts const test_shifts = boost::assign::list_of (0)(-1)(0);
	//devel::denovo_design::components::SheetList const & list = db.sheet_list( test_lengths, test_orientations, test_shifts );

	core::Real random = 0.25001;
	StructureDataOP sd_ptr = sheet.design( core::pose::Pose(), random );
	TS_ASSERT( sd_ptr );
	StructureData const & sd = *sd_ptr;

	ExtendedPoseBuilder builder;
	core::pose::PoseOP pose_ptr = builder.apply( sd );
	TS_ASSERT( pose_ptr );
	core::pose::Pose const & pose = *pose_ptr;

	sd.check_pose_consistency( pose );

	// Strands should all be in same MG
	TS_ASSERT_EQUALS( sd.segment( s1_id ).movable_group(), sd.segment( s2_id ).movable_group() );
	TS_ASSERT_EQUALS( sd.segment( s1_id ).movable_group(), sd.segment( s3_id ).movable_group() );

	// SD should have sheet idx stored
	TS_ASSERT( sd.has_data_int( sheet_id, "sheet_idx" ) );
	// sheet 55 should be selected, 223 sheets * 0.25001 = 55.75 + 1  for 1-indexing = 56
	TS_ASSERT_EQUALS( sd.get_data_int( sheet_id, "sheet_idx" ), 56 );

	// compare coordinates
	}
	*/

};

