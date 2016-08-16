// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/DeNovoArchitectTests.cxxtest.hh
/// @brief  Tests for DeNovoArchitect API
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
#include <core/import_pose/import_pose.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR("DeNovoArchitectTests");

using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::components;
/// silly class for testing
class DummyDeNovoArchitect : public DeNovoArchitect {
public:
	DummyDeNovoArchitect( std::string const & id ):
		DeNovoArchitect( id ),
		number_( 0 )
	{}

	virtual std::string
	type() const { return "Dummy"; }

	virtual DeNovoArchitectOP
	clone() const { return DeNovoArchitectOP( new DummyDeNovoArchitect( *this ) ); }

	int number_;

	virtual StructureDataOP
	design( core::pose::Pose const &, core::Real & ) const
	{
		StructureDataOP sd( new StructureData( "DummyTest" ) );
		Segment seg;
		seg.extend( "LHHHL", "XBBBX" );
		sd->add_segment( "DummySegment", seg );
		sd->set_data_int( "DummyTest", "Lebron_James_Number", number_ );
		return sd;
	}

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
	{
		number_ = tag->getOption< core::Size >( "number" );
	}

};
typedef utility::pointer::shared_ptr< DummyDeNovoArchitect > DummyDeNovoArchitectOP;

class DeNovoArchitectTests : public CxxTest::TestSuite {
	//Define Variables
public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_apply()
	{
		DummyDeNovoArchitect dummy( "ApplyTest" );
		TS_ASSERT_EQUALS( dummy.number_, 0 );

		StructureDataCOP sd = dummy.apply( core::pose::Pose() );
		TS_ASSERT( sd->has_segment( "DummySegment" ) );
		TS_ASSERT_EQUALS( sd->segment( "DummySegment" ).ss(), "LHHHL" );
		TS_ASSERT_EQUALS( sd->segment( "DummySegment" ).abego(), "XBBBX" );
		TS_ASSERT( sd->has_data_int( "DummyTest", "Lebron_James_Number" ) );
		TS_ASSERT_EQUALS( sd->get_data_int( "DummyTest", "Lebron_James_Number" ), 0 );
	}

	void test_cloning()
	{
		std::string const my_id = "MyID";
		DummyDeNovoArchitect dummy( my_id );
		dummy.number_ = 23;

		DeNovoArchitectOP base_copy = dummy.clone();
		DummyDeNovoArchitectOP copy = utility::pointer::dynamic_pointer_cast< DummyDeNovoArchitect >( base_copy );
		TS_ASSERT( copy );
		TS_ASSERT_EQUALS( copy->number_, dummy.number_ );
		TS_ASSERT_EQUALS( copy->id(), dummy.id() );
	}

	void test_parse_tag()
	{
		DummyDeNovoArchitect dummy( "UnParsedTag" );
		TS_ASSERT_EQUALS( dummy.id(), "UnParsedTag" );
		TS_ASSERT_EQUALS( dummy.number_, 0 );

		std::string const my_id = "ParsedTag";
		std::stringstream xml;
		xml << "<DummyDeNovoArchitect name=\"" << my_id << "\" number=\"" << 23 << "\" />";
		utility::tag::TagCOP tag = utility::tag::Tag::create( xml );
		basic::datacache::DataMap data;
		dummy.parse_my_tag( tag, data );

		TS_ASSERT_EQUALS( dummy.id(), my_id );
		TS_ASSERT_EQUALS( dummy.number_, 23 );
	}

};
