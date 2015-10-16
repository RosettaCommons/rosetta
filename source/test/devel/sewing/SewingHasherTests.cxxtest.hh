// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/zn_hash/ZnHash.cxxtest.hh
/// @brief  test suite for devel::zn_hash::ZnHash
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <devel/sewing/hashing/Hasher.hh>
#include <test/core/init_util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility headers
#include <utility/sql_database/DatabaseSessionManager.hh>

/// Project headers

//Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel headers
#include <devel/init.hh>

// --------------- Test Class --------------- //
using namespace devel::sewing;

class SewingHasherTests : public CxxTest::TestSuite {

public:
	
	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_key_generation() {
		//devel::sewing::SewHash hasher;
		//TS_ASSERT_EQUALS(false, hasher.hash_map().key_eq()(key2, key3));
	}

	void test_transform_coords() {
		//devel::sewing::SewHash hasher;
	}

	void test_scoring() {
	}

	void test_serialization() {
		//TS_ASSERT_EQUALS(hasher.hash_map().size(), hasher2.hash_map().size());
	}

	void test_score_one() {
        using namespace devel::sewing;

        int model_id_1=13327;
        core::Size resnum_1=18;

        int model_id_2=24789;
        core::Size resnum_2=318;

        int model_id_3=13328;
        core::Size resnum_3=18;

        std::map<int, Model> models = read_model_file("devel/sewing/inputs/test3.models");
        Model model_1 = models[model_id_1] ;
        SewResidue residue_1 = model_1.get_residue(resnum_1);

        Model model_2 = models[model_id_2] ;
        SewResidue residue_2 = model_2.get_residue(resnum_2);

        Model model_3 = models[model_id_3] ;
        SewResidue residue_3 = model_3.get_residue(resnum_3);

		Hasher hasher;
        ScoreResult result = hasher.score_one(model_1, residue_1, model_2, residue_2);
        
        ScoreResult result2 = hasher.score_one(model_3, residue_3, model_2, residue_2);

        std::cout << "Found1 " << result.second.segment_matches.size() << " matching segments" << std::endl;
        std::cout << "Found2 " << result2.second.segment_matches.size() << " matching segments" << std::endl;
        TS_ASSERT_EQUALS(result.second.segment_matches.size(), result2.second.segment_matches.size());

    }

	void test_hashing() {
        using namespace devel::sewing;
		Hasher hasher;

        utility::vector1< std::pair<core::Size, core::Size> > segments;
        segments.push_back(std::make_pair(1,14));
        segments.push_back(std::make_pair(15,28));
        segments.push_back(std::make_pair(29,42));

        core::pose::Pose bundle_1_pose;
        core::import_pose::pose_from_pdb( bundle_1_pose, "devel/sewing/inputs/bundle1.pdb" );
        Model bundle_1_model = create_model_from_pose(bundle_1_pose, segments, 1);

        core::pose::Pose bundle_2_pose;
        core::import_pose::pose_from_pdb( bundle_2_pose, "devel/sewing/inputs/bundle2.pdb" );
        Model bundle_2_model = create_model_from_pose(bundle_2_pose, segments, 2);

        hasher.insert(bundle_1_model); 
        ScoreResults results = hasher.score(bundle_2_model, 2/*segment matches*/, 10/*min segment score*/, 0/*max clash*/, true/*store atoms*/);

        TS_ASSERT_EQUALS( results.size(), 1 );

        HashResult hash_result = results.begin()->second;
        TS_ASSERT_EQUALS( hash_result.clash_count, 0 );
        TS_ASSERT_EQUALS( hash_result.segment_matches.size(), 2 );

        std::pair< SegmentPair, AtomMap > segment_match = *hash_result.segment_matches.begin();
        TS_ASSERT_EQUALS( segment_match.second.size(), 52 );//14*4 basis atoms - 4 for the aligned residue
	}

private:

};
