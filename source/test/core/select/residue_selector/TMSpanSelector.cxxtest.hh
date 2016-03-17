// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/select/residue_selector/TMSpanSelector.cxxtest.hh
/// @brief  Select residues within given transmembrane spans in a membrane protein
/// @author Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/select/residue_selector/TMSpanSelector.hh> 
#include <core/select/residue_selector/ResidueSelector.hh> 

#include <protocols/membrane/AddMembraneMover.hh> 

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh> 
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("TMSpanSelector");

class TMSpanSelectorTest : public CxxTest::TestSuite {

public:

	void setUp(){

		using namespace core::pose; 
		using namespace core::import_pose; 
		using namespace protocols::membrane; 

		core_init();

		// Setup DsbB (4 TM segments) and define transmembrane spans during add
		// membrane mover setup
		core::pose::PoseOP temp_pose = pose_from_file("core/select/residue_selector/DsbB_4tm_segs.pdb", core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb( new AddMembraneMover( "core/select/residue_selector/DsbB_4tm_segs_B.span" ) );
		add_memb->apply( *temp_pose );
		pose_ = core::pose::PoseOP( new Pose( *temp_pose ) ); 

	}

	void tearDown() {

	}

	void test_select_all_tm_segments() {

		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing selection of all transmembrane segments in the pose" );

		utility::vector1< core::Size > segments; 
		TMSpanSelectorOP tm_selector = TMSpanSelectorOP( new TMSpanSelector( true, segments ) );
		ResidueSubset subset = tm_selector->apply( *pose_ ); 

		utility::vector1< core::Size > expected_residues; 
		for ( core::Size ii = 189; ii <= 204; ++ii ) {
			expected_residues.push_back( ii ); 
		}
		for ( core::Size ii = 224; ii <= 239; ++ii ) {
			expected_residues.push_back( ii ); 
		}
		for ( core::Size ii = 246; ii <= 260; ++ii ) {
			expected_residues.push_back( ii ); 
		}
		for ( core::Size ii = 308; ii <= 322; ++ii ) {
			expected_residues.push_back( ii ); 
		}

		for ( core::Size jj = 1; jj <= expected_residues.size(); ++jj ) {
			TS_ASSERT( subset[ expected_residues[jj] ] );
		}
	}

	void test_select_some_tm_segments() {
		
		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing selection of some transmembrane segments in the pose" ); 

		TMSpanSelectorOP tm_selector( new TMSpanSelector );
		tm_selector->all( false ); 
		tm_selector->add_span( 1 ); 
		tm_selector->add_span( 2 ); 
		ResidueSubset subset = tm_selector->apply( *pose_ );

		utility::vector1< core::Size > expected_residues; 
		for ( core::Size ii = 189; ii <= 204; ++ii ) {
			expected_residues.push_back( ii ); 
		}
		for ( core::Size ii = 224; ii <= 239; ++ii ) {
			expected_residues.push_back( ii ); 
		}

		for ( core::Size jj = 1; jj <= expected_residues.size(); ++jj ) {
			TS_ASSERT( subset[ expected_residues[jj] ] );
		}
	}

	void test_select_one_tm_segment() {

		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing selection of a single transmembrane segment in the pose" ); 

		utility::vector1< core::Size > segments; 
		segments.push_back( 2 );  
		TMSpanSelectorOP tm_selector( new TMSpanSelector( false, segments ) );
		ResidueSubset subset = tm_selector->apply( *pose_ ); 

		utility::vector1< core::Size > expected_residues; 
		for ( core::Size ii = 224; ii <= 239; ++ii ) {
			expected_residues.push_back( ii ); 
		}

		for ( core::Size jj = 1; jj <= expected_residues.size(); ++jj ) {
			TS_ASSERT( subset[ expected_residues[jj] ] );
		}
	}

	void test_select_no_tm_segments() {

		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing that an error is thrown if no transmembrane segments are specified" ); 

		utility::vector1< core::Size > segments;  
		TMSpanSelectorOP tm_selector( new TMSpanSelector( false, segments ) );
		TS_ASSERT_THROWS_ANYTHING( tm_selector->apply( *pose_ ) ); 
	}

	void test_selected_out_of_span_bounds() {

		using namespace core::select::residue_selector; 
		TS_TRACE( "Testing that an error is thrown if the selected span is out of bounds" ); 

		utility::vector1< core::Size > segments;  
		segments.push_back( 5 ); 
		TMSpanSelectorOP tm_selector( new TMSpanSelector( true, segments ) );
		TS_ASSERT_THROWS_ANYTHING( tm_selector->apply( *pose_ ) ); 
	}

private: 

	core::pose::PoseCOP pose_; 

};



