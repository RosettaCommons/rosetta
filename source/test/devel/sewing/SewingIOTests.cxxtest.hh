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

class SewingIOTests : public CxxTest::TestSuite {

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

	void test_model_IO() {
        //Read file from disk and verify
        std::map<int, Model> models = read_model_file("devel/sewing/inputs/test.models");
        TS_ASSERT_EQUALS( models.size(), 2 );
        
        Model model_1 = models[1];
        TS_ASSERT_EQUALS( model_1.model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.size(), 3 );

        TS_ASSERT_EQUALS( model_1.segments_.front().model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.front().residues_.size(), 14 );

        SewResidue model_1_first_res = model_1.segments_.front().residues_.front();
        TS_ASSERT_EQUALS( model_1_first_res.basis_atoms_.size(), 4 );
        TS_ASSERT_EQUALS( model_1_first_res.residue_type_, "VAL_p:NtermProteinFull");

        SewResidue model_1_last_res = model_1.segments_.front().residues_.back();
        TS_ASSERT_EQUALS( model_1_last_res.basis_atoms_.size(), 4 );
        TS_ASSERT_EQUALS( model_1_last_res.residue_type_, "PHE");

        TS_ASSERT_EQUALS( model_1.segments_.back().model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.back().residues_.size(), 14 );

        
        //Now write and re-read (to test writing)
        write_model_file(models, "unit_test.models");
        models = read_model_file("unit_test.models");
        TS_ASSERT_EQUALS( models.size(), 2 );

        model_1 = models[1];
        TS_ASSERT_EQUALS( model_1.model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.size(), 3 );

        TS_ASSERT_EQUALS( model_1.segments_.front().model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.front().residues_.size(), 14 );

        model_1_first_res = model_1.segments_.front().residues_.front();
        TS_ASSERT_EQUALS( model_1_first_res.basis_atoms_.size(), 4 );
        TS_ASSERT_EQUALS( model_1_first_res.residue_type_, "VAL_p:NtermProteinFull");

        model_1_last_res = model_1.segments_.front().residues_.back();
        TS_ASSERT_EQUALS( model_1_last_res.basis_atoms_.size(), 4 );
        TS_ASSERT_EQUALS( model_1_last_res.residue_type_, "PHE");

        TS_ASSERT_EQUALS( model_1.segments_.back().model_id_, 1 );
        TS_ASSERT_EQUALS( model_1.segments_.back().residues_.size(), 14 );
	}

	//void test_model_writing() {
    //    Model model_1;
	//}

    //void test_model_IO() {
    //    utility::vector1<numeric::xyzVector1<core::Real>> dummy_coords;
    //    dummy_coords.push_back(numeric::xyzVector1<core::Real>(1.0, 2.2352, -13.29));
    //    dummy_coords.push_back(numeric::xyzVector1<core::Real>(-19.0, -5.55, -13.29));
    //    dummy_coords.push_back(numeric::xyzVector1<core::Real>(-4.20, 5.0, 8.81));

    //    //Model 1
    //    Model m1;
    //    m1.model_id_ = 1;
    //    for(core::Size i=1; i<=4; ++i){
    //        SewSegment s;
    //        for(core::Size j=1; j<=4; ++j){
    //            SewResidue r;
    //            for(core::Size k=1; k<=3; ++k) {
    //                SewAtom a;
    //                a.atomno = k;
    //                a.coords = dummy_coords[k];
    //                r.basis_atoms_.push_back(a);
    //            }
    //            r.num_neighbors_ = 8;
    //            s.residues_.push_back(r);
    //        }
    //        s.model_id = 1;
    //        s.segment_id = i;
    //        m1.segments_.push_back(s);
    //    }

    //    //Model 2
    //    Model m2;
    //    m2.model_id_ = 2;
    //    for(core::Size i=4; i<=1; --i){
    //        SewSegment s;
    //        for(core::Size j=5; j<=1; --j){
    //            SewResidue r;
    //            for(core::Size k=3; k<=1; --k) {
    //                SewAtom a;
    //                a.atomno = k;
    //                a.coords = dummy_coords[k];
    //                r.basis_atoms_.push_back(a);
    //            }
    //            r.num_neighbors_ = 4;
    //            s.residues_.push_back(r);
    //        }
    //        s.model_id = 2;
    //        s.segment_id = i;
    //        m2.segments_.push_back(s);
    //    }

    //    std::map<int, Model> models;
    //    models.insert(std::make_pair(1, m1));
    //    models.insert(std::make_pair(2, m2));
    //    write_model_file(models, "unit_test.models");
    //    
    //}

	void test_ascii_score_reading() {
	}

	void test_binary_score_reading() {
	}

	void test_ascii_score_writing() {
	}

	void test_binary_score_writing() {
	}

private:

};
