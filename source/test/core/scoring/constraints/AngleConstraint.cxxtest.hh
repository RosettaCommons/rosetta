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

// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <platform/types.hh>
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
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
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
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/FourPointsFunc.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
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
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/Tracer.fwd.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class AngleConstraintTests : public CxxTest::TestSuite
{

public:
	AngleConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_angle_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		HarmonicFuncOP func = new HarmonicFunc( numeric::conversions::radians( 109 ), 10 );

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );

		AngleConstraint ang_cst( at1, at2, at3, func );
		EnergyMap weights, emap;
		weights[ angle_constraint ] = 1.0;
		ang_cst.score( fourpts, weights, emap );
		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ angle_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ angle_constraint ],   0.0138748361634874, 1e-14 );

		std::cout.precision( before_precision );
	}

	void test_angle_derivatives() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->total_residue() == 2 );
		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 );
		HarmonicFuncOP func = new HarmonicFunc( numeric::conversions::radians( 109 ), 10 );
		AngleConstraint ang_cst( at1, at2, at3, func );
		EnergyMap weights, emap;
		weights[ angle_constraint ] = 1.0;
		ConformationXYZ cfunc( ubqstump->conformation() );
		ang_cst.score( cfunc, weights, emap );
		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ angle_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ angle_constraint ],   0.03488216167816781, 1e-14 );

		ScoreFunction sfxn;
		sfxn.set_weight( angle_constraint, 1.0 );
		kinematics::MoveMap movemap;
    movemap.set_bb( true );
    movemap.set_chi( true );

		for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang )
			{
				at1.rsd() = at2.rsd() = at3.rsd() = ii;
				at1.atomno() = ( rsd_type.bondangle( bondang ) ).key1();
				at2.atomno() = ( rsd_type.bondangle( bondang ) ).key2();
				at3.atomno() = ( rsd_type.bondangle( bondang ) ).key3();

				AngleConstraintOP ang_cst2 = new AngleConstraint( at1, at2, at3, func );
				ubqstump->remove_constraints();
				ubqstump->add_constraint( ang_cst2 );
				AtomDerivValidator adv( *ubqstump, sfxn, movemap );
				/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
				/// that the analytic norm matches the numeric norm to 1e-3.
				adv.simple_deriv_check( true, 1e-6 );
			}
		}


	}


};
