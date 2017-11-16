// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/components/FoldGraphTests.cxxtest.hh
/// @brief  Test suite for FoldGraph class
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Protocol Headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

// Basic/Utility Headers
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/assign.hpp>

// C++ Headers
#include <set>

static basic::Tracer TR("FoldGraphTests");


class FoldGraphTests : public CxxTest::TestSuite {
	typedef protocols::denovo_design::SegmentNames SegmentNames;
	typedef protocols::denovo_design::components::FoldGraph FoldGraph;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;
	typedef protocols::denovo_design::components::StructureDataFactory StructureDataFactory;
	//Define Variables
public:

	void setUp()
	{
		core_init();

		// Set up -run:preserver header before reading any input pdbs
		StructureDataFactory::get_instance();
	}

	void tearDown()
	{
	}

	void test_simple_fold_tree()
	{
		TR << "TestSimpleFoldTree" << std::endl;
	}

	void test_create_loops()
	{
		// Create a test StructureData
		std::string const sd_id = "TestCreateLoops";
		std::stringstream xml;
		xml << "<StructureData name=" << sd_id << " length=34 pose_length=38 >" << std::endl;
		xml << "<ResidueRange name=h1 start=2 group=1 nterm=0 cterm=1 ss=LHHHHHHHHHH abego=XAAAAAAAAAA upper_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=l1 start=12 group=3 nterm=1 cterm=1 ss=LLLL abego=BGAB lower_segment=h1 upper_segment=h2 />" << std::endl;
		xml << "<ResidueRange name=h2 start=16 group=2 nterm=1 cterm=0 ss=HHHHHHHHL abego=AAAAAAAAX lower_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=fixed start=25 nterm=0 cterm=0 group=1 ss=LEEEEELLEEEEEL abego=XBBBBBGGBBBBBX />" << std::endl;
		xml << "</StructureData>" << std::endl;
		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( xml );

		FoldGraph fg( sd );

		// list existing loops
		SegmentNames const staple_loops = boost::assign::list_of ("l1");

		// should create a loop between end of h1 and start of h2 with no cuts
		protocols::loops::LoopsOP remodel_regions = fg.create_loops( staple_loops );
		TS_ASSERT( remodel_regions );
		TS_ASSERT_EQUALS( remodel_regions->num_loop(), 1 );
		TS_ASSERT_EQUALS( remodel_regions->begin()->start(), sd.segment("l1").lower() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->stop(), sd.segment("h2").upper() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->cut(), 0 );

		// now try with h2 in same movable group as fixed
		sd.set_movable_group( "h1", 2 );
		sd.set_movable_group( "h2", 1 );

		fg = FoldGraph( sd );
		remodel_regions = fg.create_loops( staple_loops );
		TS_ASSERT( remodel_regions );
		TS_ASSERT_EQUALS( remodel_regions->num_loop(), 1 );
		TS_ASSERT_EQUALS( remodel_regions->begin()->start(), sd.segment("h1").lower() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->stop(), sd.segment("l1").upper() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->cut(), 0 );


		// now try setting a cutpoint
		sd.set_cutpoint( "l1", 2 );
		fg = FoldGraph( sd );
		remodel_regions = fg.create_loops( SegmentNames() );
		TS_ASSERT( remodel_regions );
		TS_ASSERT_EQUALS( remodel_regions->num_loop(), 1 );
		TS_ASSERT_EQUALS( remodel_regions->begin()->start(), sd.segment("l1").lower() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->stop(), sd.segment("l1").upper() );
		TS_ASSERT_EQUALS( remodel_regions->begin()->cut(), sd.segment("l1").cutpoint() );


	}

	void test_create_fold_tree()
	{
		// Create a test StructureData
		std::string const sd_id = "TestCreateLoops";
		std::stringstream xml;
		xml << "<StructureData name=" << sd_id << " length=34 pose_length=38 >" << std::endl;
		xml << "<ResidueRange name=h1 start=2 group=1 nterm=0 cterm=1 ss=LHHHHHHHHHH abego=XAAAAAAAAAA upper_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=l1 start=12 group=3 nterm=1 cterm=1 ss=LLLL abego=BGAB lower_segment=h1 upper_segment=h2 />" << std::endl;
		xml << "<ResidueRange name=h2 start=16 group=2 nterm=1 cterm=0 ss=HHHHHHHHL abego=AAAAAAAAX lower_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=fixed start=25 nterm=0 cterm=0 group=1 ss=LEEEEELLEEEEEL abego=XBBBBBGGBBBBBX />" << std::endl;
		xml << "</StructureData>" << std::endl;
		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( xml );
		FoldGraph fg( sd );

		// root fold tree at h2
		std::string root = "h2";
		SegmentNames const roots = boost::assign::list_of (root);
		core::kinematics::FoldTree const ft = fg.fold_tree( roots );

		TS_ASSERT( ft.check_fold_tree() );
		TS_ASSERT_EQUALS( ft.num_jump(), 1 );
		TS_ASSERT_EQUALS( ft.num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( ft.root(), sd.segment(root).safe() );
	}

	void test_fold_tree_with_cutpoint()
	{
		// Create a test StructureData
		std::string const sd_id = "TestCreateLoops";
		std::stringstream xml;
		xml << "<StructureData name=" << sd_id << " length=34 pose_length=38 >" << std::endl;
		xml << "<ResidueRange name=h1 start=2 group=1 nterm=0 cterm=1 ss=LHHHHHHHHHH abego=XAAAAAAAAAA upper_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=l1 start=12 group=3 nterm=1 cterm=1 cutpoint=2 ss=LLLL abego=BGAB lower_segment=h1 upper_segment=h2 />" << std::endl;
		xml << "<ResidueRange name=h2 start=16 group=2 nterm=1 cterm=0 ss=HHHHHHHHL abego=AAAAAAAAX lower_segment=l1 />" << std::endl;
		xml << "<ResidueRange name=fixed start=25 nterm=0 cterm=0 group=1 ss=LEEEEELLEEEEEL abego=XBBBBBGGBBBBBX />" << std::endl;
		xml << "</StructureData>" << std::endl;
		StructureData sd = StructureDataFactory::get_instance()->create_from_xml( xml );
		FoldGraph fg( sd );

		// root fold tree at h2
		std::string const root = "h2";
		SegmentNames const roots = boost::assign::list_of (root);
		core::kinematics::FoldTree const ft = fg.fold_tree( roots );

		// TODO: Figure out what to do with cutpoints in fold trees
		/*
		std::set< int > const wanted_cuts = boost::assign::list_of (sd.segment(root).upper())(sd.segment("l1").cutpoint());

		TS_ASSERT( ft.check_fold_tree() );
		TS_ASSERT_EQUALS( ft.num_jump(), 2 );
		TS_ASSERT_EQUALS( ft.num_cutpoint(), 2 );
		TS_ASSERT_EQUALS( ft.root(), sd.segment(root).safe() );

		TS_ASSERT_EQUALS( wanted_cuts.size(), ft.cutpoints().size() );
		for ( utility::vector1< int >::const_iterator cut=ft.cutpoints().begin(); cut!=ft.cutpoints().end(); ++cut ) {
		TS_ASSERT_DIFFERS( wanted_cuts.find( *cut ), wanted_cuts.end() );
		}
		*/
	}
};

