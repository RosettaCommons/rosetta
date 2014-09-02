// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/backrub/BackrubMover.cxxtest.hh
/// @brief  test suite for protocols::backrub::BackrubMover.cc
/// @author Colin A. Smith

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Protocol Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/backrub/BackrubMover.hh>

// Numeric Headers
#include <numeric/conversions.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <protocols/branch_angle/BranchAngleOptimizer.fwd.hh>
#include <protocols/branch_angle/BranchCoef1.fwd.hh>
#include <protocols/branch_angle/BranchCoef2.fwd.hh>
#include <protocols/branch_angle/BranchParam1.fwd.hh>
#include <protocols/branch_angle/BranchParam2.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/backrub/BackrubSegment.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key3Vector.fwd.hh>
#include <utility/keys/Key3Vector.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/IntervalSet.fwd.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>


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

		core::init::init_random_generators(1000, numeric::random::_RND_TestRun_, "mt19937");
	}

	void tearDown() {
		the_pose = 0;
	}

	void test_BackrubMover() {

		//test::UTracer UT("prococols/moves/BackrubMover.u");

		protocols::backrub::BackrubMover backrubmover;

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
