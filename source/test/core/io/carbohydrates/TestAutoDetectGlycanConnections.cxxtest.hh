// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/io/carbohydrates/TestAutoDetectGlycanConnections.cxxtest.hh
/// @brief  Test the auto-detect glycan logic of PDB IO.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>


// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
//#include <core/scoring/ScoreFunction.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

static basic::Tracer TR("TestAutoDetectGlycanConnections");


class TestAutoDetectGlycanConnections : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-include_sugars -auto_detect_glycan_connections -alternate_3_letter_codes pdb_sugar -ignore_unrecognized_res -load_PDB_components false" );

	}

	void tearDown(){

	}

	void test_empty() {
		// have to provide at least one test or cxxtest gets grumpy
		TS_ASSERT( true );
	}

	void test_auto_detect_io(){
		using namespace core::pose;
		using namespace core::import_pose;

		Pose p1;

		utility::vector1< std::string > pdb_files( {
			"1ioo",
			"1x38",
			"3uue",
			"2cl2",
			"1jnd",
			"4q56",
			"4nyq",
			"4f8x"
			});

		std::map< std::string, std::string > residue_orders ( {
			{ "1ioo", "A:1-199 A:201 A:200 B:1-199 B:201-203 B:200" },
			{ "1x38", "A:1-602 A:2211-2213 A:4981 A:4984 A:4982-4983 A:4985-4987 A:6001 A:6004 A:6002-6003 A:6005" },
			{ "3uue", "A:26-304 A:501-506" },
			{ "2cl2", "A:1-298 A:1299-1302 A:1305 A:1303-1304" },
			{ "1jnd", "A:2-141 A:161-420 A:3100-3103" },
			{ "4q56", "A:2-101 A:202-204 A:207-209 A:205-206 A:210" },
			{ "4nyq", "A:1-153 A:201-206" },
			{ "4f8x", "A:26-360 A:401-403 A:405 A:404 A:406" }
			});

		std::map< std::string, std::string > fold_trees ( {
			// Foldtree for 1ioo is a bit messed up because of the missing XYP -- The 199-200 edge is a chemical connection, instead of being a backbone one, and the 399-400 one is technically being listed as cutpoints in the ResidueTypes
			{ "1ioo", "FOLD_TREE  EDGE 1 196 -1  EDGE 28 197 -2  ND2  C1   EDGE 197 200 -1  EDGE 1 201 1  EDGE 1 202 2  EDGE 202 397 -1  EDGE 229 398 -2  ND2  C1   EDGE 398 402 -1  EDGE 400 403 -2  O6   C1   EDGE 1 404 3 " },
			{ "1x38", "FOLD_TREE  EDGE 1 602 -1  EDGE 221 603 -2  ND2  C1   EDGE 603 605 -1  EDGE 498 606 -2  ND2  C1   EDGE 606 607 -1  EDGE 606 608 -2  O4   C1   EDGE 608 611 -1  EDGE 1 612 1  EDGE 600 613 -2  ND2  C1   EDGE 613 614 -1  EDGE 613 615 -2  O4   C1   EDGE 615 616 -1  EDGE 1 617 2 " },
			{ "3uue", "FOLD_TREE  EDGE 1 279 -1  EDGE 228 280 -2  ND2  C1   EDGE 280 284 -1  EDGE 7 285 -2  OG1  C1  " },
			{ "2cl2", "FOLD_TREE  EDGE 1 298 -1  EDGE 43 299 -2  ND2  C1   EDGE 299 303 -1  EDGE 302 304 -2  O6   C1   EDGE 304 305 -1 " },
			{ "1jnd", "FOLD_TREE  EDGE 1 400 -1  EDGE 180 401 -2  ND2  C1   EDGE 401 404 -1 " },
			{ "4q56", "FOLD_TREE  EDGE 1 100 -1  EDGE 1 101 1  EDGE 1 102 2  EDGE 34 103 -2  ND2  C1   EDGE 103 106 -1  EDGE 103 107 -2  O4   C1   EDGE 107 108 -1  EDGE 1 109 3 " },
			{ "4nyq", "FOLD_TREE  EDGE 1 153 -1  EDGE 35 154 -2  ND2  C1   EDGE 154 156 -1  EDGE 79 157 -2  ND2  C1   EDGE 157 158 -1  EDGE 145 159 -2  ND2  C1  " },
			{ "4f8x", "FOLD_TREE  EDGE 1 335 -1  EDGE 311 336 -2  ND2  C1   EDGE 336 339 -1  EDGE 338 340 -2  O6   C1   EDGE 63 341 -2  ND2  C1  " }
			});


		for ( std::string const & pdb_file : pdb_files ) {
			std::string fullname( "core/io/carbohydrates/" + pdb_file + ".pdb" );
			TS_ASSERT_THROWS_NOTHING( pose_from_file( p1, fullname, PDB_file) );

			//Debug output
			//core::scoring::ScoreFunction empty;
			//empty( p1 );
			//core::io::pdb::dump_pdb( p1, fullname + ".post.pdb" );

			// Check for proper ordering
			std::string short_desc( p1.pdb_info()->short_desc() );
			TR << pdb_file << " ORDER: " << short_desc << std::endl;
			TS_ASSERT_EQUALS( short_desc, residue_orders[ pdb_file ] );

			// Check for proper sequence
			TR << pdb_file << " SEQ: " << p1.annotated_sequence() << std::endl;
			std::string sequence;
			utility::io::izstream seq_stream( "core/io/carbohydrates/" + pdb_file + ".ann_seq" );
			seq_stream >> sequence; // Should be no spaces in the sequences
			TS_ASSERT_EQUALS( p1.annotated_sequence(), sequence );

			// Check for proper FoldTree connectivity
			TR << pdb_file << " FT: '" << p1.fold_tree().to_string() << "'" << std::endl;
			TS_ASSERT_EQUALS( p1.fold_tree().to_string(), fold_trees[ pdb_file ] );
		}

	}



};



