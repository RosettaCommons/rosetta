// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/RTConstraint.cxxtest.hh
/// @brief  test suite for rt constraint
/// @author atom-moyer (apmoyer@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.hh>

#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/constraints/RTConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>


//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/FourPointsFunc.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/Tracer.fwd.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.RTConstraint.cxxtest");

using namespace core;

class RTConstraintTests : public CxxTest::TestSuite
{

public:
	RTConstraintTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_rt_constraint_zero_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();

		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 1.0 );

		TS_ASSERT( ubqstump->size() == 2 );

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );

		AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
		AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );

		StubID stub_id1( atom_id1, atom_id2, atom_id3 );
		StubID stub_id2( atom_id4, atom_id5, atom_id6 );

		kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz(atom_id1),
			ubqstump->xyz(atom_id2),
			ubqstump->xyz(atom_id3)),
			kinematics::Stub(ubqstump->xyz(atom_id4),
			ubqstump->xyz(atom_id5),
			ubqstump->xyz(atom_id6)));

		RTConstraintOP rt_cst( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );

		ubqstump->add_constraint( rt_cst );

		sfxn.score( *ubqstump );

		Real atom_pair_score( ubqstump->energies().total_energies().get( atom_pair_constraint ) );

		TS_ASSERT_DELTA( atom_pair_score, 0.0, 1e-3 );

	}

	void test_rt_constraint_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();

		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 1.0 );

		TS_ASSERT( ubqstump->size() == 2 );

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.866025 ) );

		AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
		AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );

		StubID stub_id1( atom_id1, atom_id2, atom_id3 );
		StubID stub_id2( atom_id4, atom_id5, atom_id6 );

		kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz( atom_id1 ),
			ubqstump->xyz( atom_id2 ),
			ubqstump->xyz( atom_id3 )),
			kinematics::Stub(ubqstump->xyz( atom_id4 ) - 0.5,
			ubqstump->xyz( atom_id5 ) - 0.5,
			ubqstump->xyz( atom_id6 ) - 0.5));

		RTConstraintOP rt_cst( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );

		ubqstump->add_constraint( rt_cst );

		sfxn.score( *ubqstump );

		Real atom_pair_score( ubqstump->energies().total_energies().get( atom_pair_constraint ) );

		TS_ASSERT_DELTA( atom_pair_score, 1.0, 1e-3 );

	}

	// void test_rt_constraint_zero_derivative() {
	//  using namespace core;
	//  using namespace core::id;
	//  using namespace core::scoring;
	//  using namespace core::scoring::constraints;
	//
	//  core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
	//
	//  ScoreFunction sfxn;
	//  sfxn.set_weight( atom_pair_constraint, 1.0 );
	//
	//  kinematics::MoveMap movemap;
	//  movemap.set_bb( true );
	//
	//  TS_ASSERT( ubqstump->size() == 2 );
	//
	//  core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );
	//
	//  AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
	//  AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );
	//
	//  StubID stub_id1( atom_id1, atom_id2, atom_id3 );
	//  StubID stub_id2( atom_id4, atom_id5, atom_id6 );
	//
	//  kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz( atom_id1 ),
	//                                             ubqstump->xyz( atom_id2 ),
	//                        ubqstump->xyz( atom_id3 )),
	//               kinematics::Stub(ubqstump->xyz( atom_id4 ),
	//                                  ubqstump->xyz( atom_id5 ),
	//                                             ubqstump->xyz( atom_id6 )));
	//
	//  RTConstraintOP rt_cst( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );
	//
	//  ubqstump->add_constraint( rt_cst );
	//
	//  AtomDerivValidator adv( *ubqstump, sfxn, movemap );
	//
	//  adv.simple_deriv_check( true, 1e-6 );
	// }
	//
	// void test_rt_constraint_derivative() {
	//  using namespace core;
	//  using namespace core::id;
	//  using namespace core::scoring;
	//  using namespace core::scoring::constraints;
	//
	//  core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
	//
	//  ScoreFunction sfxn;
	//  sfxn.set_weight( atom_pair_constraint, 1.0 );
	//
	//  kinematics::MoveMap movemap;
	//  movemap.set_bb( true );
	//
	//  TS_ASSERT( ubqstump->size() == 2 );
	//
	//  core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );
	//
	//  AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
	//  AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );
	//
	//  StubID stub_id1( atom_id1, atom_id2, atom_id3 );
	//  StubID stub_id2( atom_id4, atom_id5, atom_id6 );
	//
	//  kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz( atom_id1 ),
	//                       ubqstump->xyz( atom_id2 ),
	//                       ubqstump->xyz( atom_id3 )),
	//               kinematics::Stub(ubqstump->xyz( atom_id4 ) - 0.5,
	//                       ubqstump->xyz( atom_id5 ) - 0.5,
	//                       ubqstump->xyz( atom_id6 ) - 0.5));
	//
	//  RTConstraintOP rt_cst( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );
	//
	//  ubqstump->add_constraint( rt_cst );
	//
	//  AtomDerivValidator adv( *ubqstump, sfxn, movemap );
	//
	//  adv.simple_deriv_check( true, 1e-6 );
	// }


	void test_rt_constraint_clone() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );

		AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
		AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );

		StubID stub_id1( atom_id1, atom_id2, atom_id3 );
		StubID stub_id2( atom_id4, atom_id5, atom_id6 );

		kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz( atom_id1 ),
			ubqstump->xyz( atom_id2 ),
			ubqstump->xyz( atom_id3 )),
			kinematics::Stub(ubqstump->xyz( atom_id4 ),
			ubqstump->xyz( atom_id5 ),
			ubqstump->xyz( atom_id6 )));

		RTConstraintOP rt_cst( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );
		ConstraintOP cloned_cst = rt_cst->clone();
		RTConstraintOP cloned_rtc = utility::pointer::dynamic_pointer_cast< RTConstraint > ( cloned_cst );

		// ensure the dynamic cast succeeds
		TS_ASSERT( cloned_rtc );

		// Make sure that the clone isn't the same as the original -- of course, right?
		TS_ASSERT_DIFFERS( rt_cst, cloned_cst );
		TS_ASSERT_DIFFERS( rt_cst, cloned_rtc );

		// check mutual equality; a == b and b == a
		TS_ASSERT( *rt_cst == *cloned_cst );
		TS_ASSERT( *cloned_cst == *rt_cst );

		// clone() should perform a deep copy of the internal func object, verifiable by looking
		// at the func pointers and making sure they point at different objects.
		TS_ASSERT_DIFFERS( & rt_cst->get_func(), & cloned_cst->get_func() );
	}

	void test_rt_constraint_equality_operator() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();

		TS_ASSERT( ubqstump->size() == 2 );

		AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
		AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );

		StubID stub_id1( atom_id1, atom_id2, atom_id3 );
		StubID stub_id2( atom_id4, atom_id5, atom_id6 );

		core::scoring::func::HarmonicFuncOP func1( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );

		core::scoring::func::HarmonicFuncOP func2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );

		kinematics::RT rt_target1(kinematics::Stub(ubqstump->xyz( atom_id1 ),
			ubqstump->xyz( atom_id2 ),
			ubqstump->xyz( atom_id3 )),
			kinematics::Stub(ubqstump->xyz( atom_id4 ),
			ubqstump->xyz( atom_id5 ),
			ubqstump->xyz( atom_id6 )));

		kinematics::RT rt_target2(kinematics::Stub(ubqstump->xyz( atom_id4 ),
			ubqstump->xyz( atom_id5 ),
			ubqstump->xyz( atom_id6 )),
			kinematics::Stub(ubqstump->xyz( atom_id1 ),
			ubqstump->xyz( atom_id2 ),
			ubqstump->xyz( atom_id3 )));

		RTConstraintOP rt_cst1( new RTConstraint( stub_id1, stub_id2, rt_target1, func1 ) );

		RTConstraintOP rt_cst2( new RTConstraint( stub_id2, stub_id1, rt_target1, func1 ) );

		RTConstraintOP rt_cst3( new RTConstraint( stub_id1, stub_id2, rt_target2, func1 ) );

		RTConstraintOP rt_cst4( new RTConstraint( stub_id1, stub_id2, rt_target1, func2 ) );

		TS_ASSERT( *rt_cst1 == *rt_cst1 )
			TS_ASSERT( *rt_cst2 == *rt_cst2 )
			TS_ASSERT( *rt_cst3 == *rt_cst3 )
			TS_ASSERT( *rt_cst4 == *rt_cst4 )

			TS_ASSERT( *rt_cst1 != *rt_cst2 );
		TS_ASSERT( *rt_cst2 != *rt_cst1 );

		TS_ASSERT( *rt_cst1 != *rt_cst3 );
		TS_ASSERT( *rt_cst3 != *rt_cst1 );

		TS_ASSERT( *rt_cst1 != *rt_cst4 );
		TS_ASSERT( *rt_cst4 != *rt_cst1 );

		TS_ASSERT( *rt_cst2 != *rt_cst3 );
		TS_ASSERT( *rt_cst3 != *rt_cst2 );

		TS_ASSERT( *rt_cst2 != *rt_cst4 );
		TS_ASSERT( *rt_cst4 != *rt_cst2 );

		TS_ASSERT( *rt_cst3 != *rt_cst4 );
		TS_ASSERT( *rt_cst4 != *rt_cst3 );
	}


	void test_rt_constraint_serialize() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 0.0, 0.5 ) );

		AtomID atom_id1( 1, 1 ); AtomID atom_id2( 2, 1 ); AtomID atom_id3( 3, 1 );
		AtomID atom_id4( 1, 2 ); AtomID atom_id5( 2, 2 ); AtomID atom_id6( 3, 2 );

		StubID stub_id1( atom_id1, atom_id2, atom_id3 );
		StubID stub_id2( atom_id4, atom_id5, atom_id6 );

		kinematics::RT rt_target(kinematics::Stub(ubqstump->xyz( atom_id1 ),
			ubqstump->xyz( atom_id2 ),
			ubqstump->xyz( atom_id3 )),
			kinematics::Stub(ubqstump->xyz( atom_id4 ),
			ubqstump->xyz( atom_id5 ),
			ubqstump->xyz( atom_id6 )));

		RTConstraintOP instance( new RTConstraint( stub_id1, stub_id2, rt_target, func ) );
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		RTConstraintOP instance2( new RTConstraint() );
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
