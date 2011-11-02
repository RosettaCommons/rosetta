// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/BackrubMover.cxxtest.hh
/// @brief  test suite for protocols::moves::BackrubMover.cc
/// @author Colin A. Smith

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Core Headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Protocol Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/moves/BackrubMover.hh>

// Numeric Headers
#include <numeric/conversions.hh>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <utility/keys/Key3Vector.hh>


using namespace core;
using namespace core::pose;
using namespace protocols::moves;

static basic::Tracer TR("protocols.moves.BackrubMover.cxxtest");

class BackrubMoverTest : public CxxTest::TestSuite {

public:

	chemical::ResidueTypeSetCAP residue_set;
	pose::PoseOP the_pose;

	BackrubMoverTest() {}

	void setUp() {
		core_init();
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

		the_pose = new Pose;
		core::import_pose::pose_from_pdb( *the_pose, "protocols/moves/test_in.pdb" );

		core::init_random_generators(1000, numeric::random::_RND_TestRun_, "mt19937");
	}

	void tearDown() {
		the_pose = 0;
	}

	void test_BackrubMover() {

		//test::UTracer UT("prococols/moves/BackrubMover.u");

		BackrubMover backrubmover;

		backrubmover.set_input_pose(the_pose);

		// These segment additions should succeed
		TR << "Attempting 3 residue segment" << std::endl;
		TS_ASSERT( backrubmover.add_segment(id::AtomID(the_pose->residue_type(2).atom_index(" CA "), 2),
		                                    id::AtomID(the_pose->residue_type(2).atom_index(" CA "), 4)) );

		TR << "Attempting 3 residue segment (reversed)" << std::endl;
		TS_ASSERT( backrubmover.add_segment(id::AtomID(the_pose->residue_type(3).atom_index(" CA "), 3),
		                                    id::AtomID(the_pose->residue_type(1).atom_index(" CA "), 1)) );

		TR << "Attempting 2 residue segment" << std::endl;
		TS_ASSERT( backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" CA "), 1),
		                                    id::AtomID(the_pose->residue_type(2).atom_index(" CA "), 2)) );

		TR << "Attempting 3 atom segment" << std::endl;
		TS_ASSERT( backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" CA "), 1),
		                                    id::AtomID(the_pose->residue_type(2).atom_index(" N  "), 2)) );

		TR << "Attempting 3 residue segment (duplicate)" << std::endl;
		TS_ASSERT( backrubmover.add_segment(id::AtomID(the_pose->residue_type(2).atom_index(" CA "), 2),
		                                    id::AtomID(the_pose->residue_type(4).atom_index(" CA "), 4)) );

		// These segment additions should fail
		TR << "Attempting 3 residue segment (invalid)" << std::endl;
		TS_ASSERT(!backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" OD1"), 1),
		                                    id::AtomID(the_pose->residue_type(3).atom_index(" CA "), 3)) );

		TR << "Attempting 2 atom segment (short)" << std::endl;
		TS_ASSERT(!backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" CA "), 1),
		                                    id::AtomID(the_pose->residue_type(1).atom_index(" C  "), 1)) );

		TR << "Attempting 2 atom segment (short, reversed)" << std::endl;
		TS_ASSERT(!backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" C  "), 1),
		                                    id::AtomID(the_pose->residue_type(1).atom_index(" CA "), 1)) );

		TR << "Attempting 1 atom segment (ultra short)" << std::endl;
		TS_ASSERT(!backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" C  "), 1),
		                                    id::AtomID(the_pose->residue_type(1).atom_index(" C  "), 1)) );

		TR << "Attempting 2 residue segment (N atoms)" << std::endl;
		TS_ASSERT(!backrubmover.add_segment(id::AtomID(the_pose->residue_type(1).atom_index(" N  "), 1),
		                                    id::AtomID(the_pose->residue_type(2).atom_index(" N  "), 2)) );

		// This section needs actual asserts/UTracer output
		int segid;

		utility::vector0<Real> start_constants;
		utility::vector0<Real> end_constants;

		TR << "Rotating segment 1" << std::endl;
		segid = 1;
		backrubmover.rotate_segment(*the_pose, segid, numeric::conversions::radians(20.), start_constants, end_constants);

		// This doesn't work with BranchAngleOptimizer
		//TR << "Rotating proline CB-CD segment" << std::endl;
		//start_constants.clear();
		//end_constants.clear();
		//segid = backrubmover.add_segment(id::AtomID(the_pose->residue_type(19).atom_index(" CB "), 19),
		//																 id::AtomID(the_pose->residue_type(19).atom_index(" CD "), 19));
		//backrubmover.rotate_segment(*the_pose, segid, numeric::conversions::radians(20.), start_constants, end_constants);

		TR << "Rotating lysine side chain segment" << std::endl;
		start_constants.clear();
		end_constants.clear();
		segid = backrubmover.add_segment(id::AtomID(the_pose->residue_type(42).atom_index(" CB "), 42),
																		 id::AtomID(the_pose->residue_type(42).atom_index("2HZ "), 42));
		backrubmover.rotate_segment(*the_pose, segid, numeric::conversions::radians(40.), start_constants, end_constants);

		//backrubmover.branchopt().read_database();
		backrubmover.apply(*the_pose);
		//backrubmover.branchopt().write_database();

		//the_pose->dump_pdb("test.pdb");
	}
};
