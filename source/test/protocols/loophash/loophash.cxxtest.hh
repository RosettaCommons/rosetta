// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/loophash.cxxtest.hh
/// @brief test suite for protocols/loophash/loophash
/// @author Mike Tyka (mike.tyka@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <basic/Tracer.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

static basic::Tracer TR("protocols.loophash.loophash.cxxtest");

namespace {

using namespace core::sequence;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::pose::make_pose_from_sequence;
using namespace protocols;
using namespace protocols::loophash;

class LoophashTest: public CxxTest::TestSuite {
 public:

	void setUp() {
		core_init();
	}
	
	
	void test_database() {
		TR << "Testing..." << std::endl;
		core::pose::Pose	sample_pose;

		std::string seq = "TAKESMEFCKANDSMSHITMAKEAFCKSHITSTACKAFCKSHITSTACK"; 
		make_pose_from_sequence(sample_pose, seq ,"fa_standard");

		// make an alpha helix
		for( core::Size ir = 1; ir < sample_pose.total_residue(); ir ++ ) {
			sample_pose.set_phi( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
			sample_pose.set_psi( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
			sample_pose.set_omega( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
		}
		
		TR << "Importing test pose ..." << std::endl;

		utility::vector1 < core::Size > loop_sizes;
		loop_sizes.push_back( 10 );
		loophash::LoopHashLibraryOP loop_hash_library( new loophash::LoopHashLibrary( loop_sizes, 0 , 1 ) );
		
		TS_ASSERT(loop_hash_library->test_saving_library( sample_pose, 14, true ));

		TR << "Save Library: " << std::endl;

		loop_hash_library->save_db();

		TR << "Reloading Library: " << std::endl;

		loophash::LoopHashLibraryOP read_loop_hash_library( new loophash::LoopHashLibrary( loop_sizes, 0 , 1 ) );
		read_loop_hash_library->load_db();
		
		TR << "Retesting Library: " << std::endl;

		TS_ASSERT(read_loop_hash_library->test_saving_library( sample_pose, 14, false  ) );
		
		TR << "Done test" << std::endl;

}

};

}  // anonymous namespace
