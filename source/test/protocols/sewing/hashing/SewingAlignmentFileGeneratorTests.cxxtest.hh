// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/sewing/SewingAlignmentFileGeneratorTests.cxxtest.hh
/// @brief  a unit test to ensure proper handling of the BasisMap
/// @author Minnie Langlois (minnie@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/sewing/hashing/AlignmentFileGenerator.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

using namespace protocols::sewing;

static basic::Tracer TR("SewingAlignmentFileGeneratorTests");


class SewingAlignmentFileGeneratorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_load_edge_and_model_file(){
		EdgeFileReaderOP edge_file_reader = EdgeFileReaderOP( new EdgeFileReader(test_edge_file_name) );
		try {
			AlignmentFileGenerator alignment_file_generator = AlignmentFileGenerator(edge_file_reader, test_model_file_name);
			TS_ASSERT_EQUALS( alignment_file_generator.version(), edge_file_reader->version() );
			Model model_to_test = (*alignment_file_generator.model_map())[ test_model_id ];
			TR << model_to_test << std::endl;
		} catch (...) {
			TS_FAIL("AFG could not import files");
		}
	}

	void test_generating_alignments(){

		EdgeFileReaderOP edge_file_reader = EdgeFileReaderOP( new EdgeFileReader(test_edge_file_name) );
		try {
			AlignmentFileGenerator alignment_file_generator = AlignmentFileGenerator(edge_file_reader, test_model_file_name);
			TS_ASSERT_EQUALS( alignment_file_generator.version(), edge_file_reader->version() );
			Model model_to_test = (*alignment_file_generator.model_map())[ test_model_id ];
			core::Size segment_id = 3;
			SewSegment segment_to_test = model_to_test.find_segment(segment_id);
			TR << "Time to Calculate Alignments" << std::endl;
			BasisMap::iterator alignments = alignment_file_generator.get_alignments(segment_to_test);
			std::pair< int , core::Size > test_seg1 = std::make_pair( 101555, 25 );
			std::pair< int , core::Size > test_seg2 = std::make_pair( 121140, 3 );
			std::pair< int , core::Size > test_seg3 = std::make_pair( 130698, 9 );
			TS_ASSERT( alignments->second.find( test_seg1 ) != alignments->second.end() );
			TS_ASSERT( alignments->second.find( test_seg2 ) != alignments->second.end() );
			TS_ASSERT( alignments->second.find( test_seg3 ) != alignments->second.end() );

			std::pair< int, core::Size > segment_key = std::make_pair( test_model_id, segment_id );
			utility::pointer::shared_ptr< BasisMap > basis_map = alignment_file_generator.basis_map();
			TS_ASSERT( basis_map->find( segment_key ) != basis_map->end() );
		} catch (...) {
			TS_FAIL("AFG could not import files");
		}

	}




private:
	std::string test_edge_file_name = "protocols/sewing/inputs/test.edge";
	std::string test_model_file_name = "protocols/sewing/inputs/5resmaxL_18resminH_helical_smotifs_20150924_SLG_wiggins.models";
	int test_model_id = 131885;

};



