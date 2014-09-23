// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/AngleConstraint.cxxtest.hh
/// @brief  test suite for angle constraints
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <numeric/conversions.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
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


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.AtomPairConstraint.cxxtest");

using namespace core;

class AtomPairConstraintTests : public CxxTest::TestSuite
{

public:
	AtomPairConstraintTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_atom_pair_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::scoring::func::FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );

		AtomID at1( 1, 1), at2( 2, 1 );

		AtomPairConstraint atom_pair_cst( at1, at2, func );
		EnergyMap weights, emap;
		weights[ atom_pair_constraint ] = 1.0;
		atom_pair_cst.score( fourpts, weights, emap );
		//Size before_precision = std::cout.precision();
		//std::cout.precision( 16 );
		//std::cout << "Atom pair constraint func: " << emap[ atom_pair_constraint ] << std::endl;
		//std::cout.precision( before_precision );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ],   0.1599999999999999, 1e-12 );

	}

	void test_atom_pair_derivatives() {

		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->total_residue() == 2 );

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );
		AtomID at1, at2;
		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 1.0 );
		kinematics::MoveMap movemap;
    movemap.set_bb( true );
    movemap.set_chi( true );

		for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size jj = 1; jj <= rsd_type.natoms(); ++jj ) {
				for ( Size kk = jj + 1; kk <= rsd_type.natoms(); ++kk ) {
					if ( rsd_type.path_distance( jj, kk ) != 1 ) continue;

					at1.rsd() = at2.rsd() = ii;
					at1.atomno() = jj;
					at2.atomno() = kk;

					AtomPairConstraintOP atom_pair_cst( new AtomPairConstraint( at1, at2, func ) );
					ubqstump->remove_constraints();
					ubqstump->add_constraint( atom_pair_cst );

					AtomDerivValidator adv( *ubqstump, sfxn, movemap );
					/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
					/// that the analytic norm matches the numeric norm to 1e-3.
					adv.simple_deriv_check( true, 1e-6 );
				}
			}
		}
	}


};
