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

// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/FourPointsFunc.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.hh>

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
#include <core/scoring/constraints/FourPointsFunc.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
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

		FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		HarmonicFuncOP func = new HarmonicFunc( 1.2, 0.5 );

		AtomID at1( 1, 1), at2( 2, 1 );

		AtomPairConstraint atom_pair_cst( at1, at2, func );
		EnergyMap weights, emap;
		weights[ atom_pair_constraint ] = 1.0;
		atom_pair_cst.score( fourpts, weights, emap );
		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		//std::cout << "Atom pair constraint func: " << emap[ atom_pair_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ],   0.1599999999999999, 1e-12 );

		std::cout.precision( before_precision );
	}

	void test_atom_pair_derivatives() {

		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->total_residue() == 2 );
		ConformationXYZ cfunc( ubqstump->conformation() );
		HarmonicFuncOP func = new HarmonicFunc( 1.2, 0.5 );
		EnergyMap weights, emap;
		weights[ atom_pair_constraint ] = 0.5;

		AtomID at1, at2;

		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );

		/*for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size jj = 1; jj <= rsd_type.natoms(); ++jj ) {
				for ( Size kk = jj + 1; kk <= rsd_type.natoms(); ++kk ) {
					if ( rsd_type.path_distance( jj, kk ) != 1 ) continue;

					at1.rsd() = at2.rsd() = ii;
					at1.atomno() = jj;
					at2.atomno() = kk;

					std::cout <<"{\n"
					"at1.rsd() = at2.rsd() = " << ii << ";\n"
					"at1.atomno() = " << jj << ";\n"
					"at2.atomno() = " << kk << ";\n"
					"AtomPairConstraint atom_pair_cst( at1, at2, func );\n"
					"atom_pair_cst.score( cfunc, weights, emap );\n";

					AtomPairConstraint atom_pair_cst( at1, at2, func );
					atom_pair_cst.score( cfunc, weights, emap );
					std::cout << "TS_ASSERT_DELTA( emap[ atom_pair_constraint ], " << emap[ atom_pair_constraint ] << ", 1e-12 );\n";

					{ /// Scope atom 1
					Vector f1( 0.0 ), f2( 0.0 );

					std::cout << "{\n"
					"//Atom 1\n"
					"Vector f1( 0.0 ), f2( 0.0 );\n"
					"atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

					atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

					std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
					std::cout << "}\n";
					}/// Scope atom 1

					{/// Scope atom 2
					Vector f1( 0.0 ), f2( 0.0 );

					std::cout << "{\n"
					"//Atom 2\n"
					"Vector f1( 0.0 ), f2( 0.0 );\n"
					"atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

					atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

					std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
					std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
					std::cout << "}\n";
					}/// Scope atom 2

					std::cout  << "}\n";
				}
			}
		}*/

		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 1;
		at2.atomno() = 2;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 0.2996150466220956, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 2.228935171593635, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -6.716164828814115, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 39.45555439042037, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.7978316702510624, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.7302314076878896, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1693720864220131, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -2.228935171593635, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 6.716164828814115, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -39.45555439042037, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.7978316702510624, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.7302314076878896, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1693720864220131, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 1;
		at2.atomno() = 9;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 0.4596150466220961, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -9.971147683700133, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 9.153384040527197, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 18.74292485167656, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.04652283824034392, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.7271206214297037, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3303501587917678, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 9.971147683700133, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -9.153384040527197, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -18.74292485167656, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.04652283824034392, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.7271206214297037, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3303501587917678, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 1;
		at2.atomno() = 10;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 0.6196150466220969, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 17.54642936168729, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -20.78576254368344, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 10.74093350943225, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.1722750327284535, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.2389266444724264, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.7437979373859276, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -17.54642936168729, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 20.78576254368344, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -10.74093350943225, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.1722750327284535, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.2389266444724264, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.7437979373859276, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 1;
		at2.atomno() = 11;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 0.7796150466220998, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -5.132037470885587, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 4.270463710359765, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 13.76529304128641, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.7487898879491324, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.1656051178241019, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.2277908001996636, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 5.132037470885587, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -4.270463710359765, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -13.76529304128641, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.7487898879491324, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.1656051178241019, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.2277908001996636, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 2;
		at2.atomno() = 3;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 1.26412277091791, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 12.61279174275054, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -14.621119426229, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 14.17273753120041, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.5818411457865235, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -1.10253051736364, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.6196113592687997, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -12.61279174275054, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 14.621119426229, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -14.17273753120041, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5818411457865235, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 1.10253051736364, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.6196113592687997, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 2;
		at2.atomno() = 5;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 1.637868990258366, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 17.88403308161734, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -19.87621872441271, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 12.44628202735369, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.9371158342763095, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.4328273307359395, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.6553314369679218, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -17.88403308161734, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 19.87621872441271, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -12.44628202735369, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.9371158342763095, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.4328273307359395, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.6553314369679218, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 2;
		at2.atomno() = 12;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 1.686423994162367, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 10.18477161994137, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -9.663049339064449, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -7.722075480554182, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.16593294927596, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.1334509038530506, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3858459902880807, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -10.18477161994137, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 9.663049339064449, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 7.722075480554182, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.16593294927596, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.1334509038530506, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3858459902880807, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 3;
		at2.atomno() = 4;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 1.690085615096335, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 1.979362224789496, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.59998575371082, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -3.015733521001755, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.09571591465004672, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01731346452046052, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.07200827289191598, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -1.979362224789496, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.59998575371082, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 3.015733521001755, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.09571591465004672, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01731346452046052, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.07200827289191598, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 5;
		at2.atomno() = 6;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 2.061117005120244, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 29.9749967642661, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -29.48280390660736, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -5.261703904592808, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.1951382266504076, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01619404370542746, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.202407745128015, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -29.9749967642661, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 29.48280390660736, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 5.261703904592808, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.1951382266504076, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01619404370542746, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 1.202407745128015, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 5;
		at2.atomno() = 13;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 2.110374368720243, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 3.187165959091696, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.9916364200107523, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -15.17243009998439, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.3504599108867839, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.2569682827780024, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.09041353276667052, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -3.187165959091696, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.9916364200107523, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 15.17243009998439, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.3504599108867839, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.2569682827780024, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.09041353276667052, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 5;
		at2.atomno() = 14;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 2.158112248820244, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 2.022054938235309, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -3.25722247174852, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 8.293190322866021, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.07689858741359255, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.4064362526965675, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1408818659294649, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -2.022054938235309, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 3.25722247174852, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -8.293190322866021, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.07689858741359255, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.4064362526965675, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1408818659294649, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 6;
		at2.atomno() = 7;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 3.858142414191437, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 33.46771850990217, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -37.7758937980328, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 17.64738068514713, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 2.003731726938946, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 1.26870153617146, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.084239936572723, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -33.46771850990217, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 37.7758937980328, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -17.64738068514713, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -2.003731726938946, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -1.26870153617146, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 1.084239936572723, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 6;
		at2.atomno() = 15;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 3.907012590547435, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -2.903348863584868, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08630748258078799, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 14.75539691379534, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.384705957498443, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.2047728793588161, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.07449899038442087, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 2.903348863584868, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08630748258078799, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -14.75539691379534, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.384705957498443, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.2047728793588161, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.07449899038442087, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 6;
		at2.atomno() = 16;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 3.955581698003436, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -1.349121048353339, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 3.427032493267289, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -9.932170207191731, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02576304864602685, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.4170172995910519, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1403896968806841, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 1.349121048353339, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -3.427032493267289, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 9.932170207191731, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02576304864602685, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.4170172995910519, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1403896968806841, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 7;
		at2.atomno() = 8;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.359163839320615, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 54.15600752841708, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -50.2501957484232, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -15.58414976601443, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.6834600844590122, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.03304932710150248, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -2.268505812246937, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -54.15600752841708, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 50.2501957484232, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 15.58414976601443, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.6834600844590122, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.03304932710150248, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 2.268505812246937, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 8;
		at2.atomno() = 17;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.407323897436615, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -7.467995891948426, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 8.32611662383488, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -2.247176579474113, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.3007479781562764, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.2031317923927685, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.2468367092192928, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 7.467995891948426, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -8.32611662383488, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 2.247176579474113, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.3007479781562764, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.2031317923927685, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.2468367092192928, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 8;
		at2.atomno() = 18;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.454979660640616, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -2.513174809027404, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.952602071637024, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 14.20876543863983, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.3834097623079768, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.2050585224954112, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.03963596012393156, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 2.513174809027404, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.952602071637024, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -14.20876543863983, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.3834097623079768, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.2050585224954112, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.03963596012393156, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 1;
		at1.atomno() = 8;
		at2.atomno() = 19;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.502877609376617, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04649499590469235, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 2.947603408002908, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -9.128432673577406, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.04206915899285039, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.4146692593308754, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.133683903936206, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04649499590469235, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -2.947603408002908, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 9.128432673577406, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.04206915899285039, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.4146692593308754, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.133683903936206, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 1;
		at2.atomno() = 2;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.854775296475186, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 10.85820029293795, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -12.03105038969595, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 14.7794857603853, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.4082613713868331, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.9917183992328702, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.5073539372574232, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -10.85820029293795, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 12.03105038969595, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -14.7794857603853, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.4082613713868331, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.9917183992328702, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5073539372574232, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 1;
		at2.atomno() = 11;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 5.999175296475183, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 13.33544653938662, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -10.58385262607794, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -17.57992547224129, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.5915963872777007, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.04371725071348294, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.4750815893612565, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -13.33544653938662, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 10.58385262607794, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 17.57992547224129, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5915963872777007, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.04371725071348294, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.4750815893612565, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 2;
		at2.atomno() = 3;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 6.414127960041218, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 31.26628800194963, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -32.11005262001231, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 23.69574249179839, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.6348218828341818, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.1963715690900397, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.103743646954364, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -31.26628800194963, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 32.11005262001231, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -23.69574249179839, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.6348218828341818, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.1963715690900397, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 1.103743646954364, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 2;
		at2.atomno() = 6;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 6.79016128082342, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -27.03482651601706, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 21.33258214785209, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 27.39667148338699, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.09524206804440978, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.9174171853508454, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.8083365262230617, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 27.03482651601706, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -21.33258214785209, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -27.39667148338699, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.09524206804440978, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.9174171853508454, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.8083365262230617, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 2;
		at2.atomno() = 12;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 6.838675748423418, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -3.631082545286514, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.514473709826196, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 13.73602483531982, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.4207447923236222, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.05681900250636897, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1174874082049786, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 3.631082545286514, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.514473709826196, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -13.73602483531982, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.4207447923236222, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.05681900250636897, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1174874082049786, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 3;
		at2.atomno() = 4;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 6.854329226092011, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.9782009856112719, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.935390925589862, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 5.975575744295679, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.2447661757538175, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.04538579291305577, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02536847813480858, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.9782009856112719, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.935390925589862, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -5.975575744295679, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.2447661757538175, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.04538579291305577, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02536847813480858, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 3;
		at2.atomno() = 5;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 6.854329226092011, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 3.843160088072e-14, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -3.26635951127313e-14, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -9.141990610162079e-15, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -8.640048923386789e-16, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -6.181128162996661e-16, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.423683852976134e-15, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -3.843160088072e-14, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 3.26635951127313e-14, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 9.141990610162079e-15, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 8.640048923386789e-16, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 6.181128162996661e-16, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 1.423683852976134e-15, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 6;
		at2.atomno() = 7;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.198960131573291, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 8.75229814991321, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -10.25159295190116, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 25.84848119527924, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.1171333276510698, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -1.099009342659024, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3962093767526095, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -8.75229814991321, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 10.25159295190116, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -25.84848119527924, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.1171333276510698, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 1.099009342659024, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3962093767526095, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 6;
		at2.atomno() = 13;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.24727832616929, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 9.002869817473853, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -9.054256293894328, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 11.11669532454315, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.3093961957210594, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.06692173777520831, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3050708990881926, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -9.002869817473853, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 9.054256293894328, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -11.11669532454315, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.3093961957210594, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.06692173777520831, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3050708990881926, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 6;
		at2.atomno() = 14;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.295769889433291, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 5.481160810058528, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -3.792549375466275, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -11.08092714758593, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.3981255079579295, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.03448025535217566, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1851308860241672, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -5.481160810058528, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 3.792549375466275, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 11.08092714758593, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.3981255079579295, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.03448025535217566, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1851308860241672, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 7;
		at2.atomno() = 8;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.729553573609747, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -33.93170000616485, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 26.09090293595613, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 26.13327531477016, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.08268789620991242, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.8751135682215708, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.9810574352405224, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 33.93170000616485, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -26.09090293595613, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -26.13327531477016, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.08268789620991242, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.8751135682215708, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.9810574352405224, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 7;
		at2.atomno() = 15;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.778211757005747, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -10.40082496754611, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 10.08003001982737, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -11.26126436723378, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.2255056941265826, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.1542839721864666, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3463760549270832, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 10.40082496754611, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -10.08003001982737, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 11.26126436723378, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.2255056941265826, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.1542839721864666, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3463760549270832, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 7;
		at2.atomno() = 16;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.826022203341747, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -3.638451200252674, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.529001987881173, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 14.54237267688264, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.4205252768412818, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.04748464748335542, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.110206556678562, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 3.638451200252674, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.529001987881173, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -14.54237267688264, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.4205252768412818, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.04748464748335542, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.110206556678562, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 8;
		at2.atomno() = 9;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.829497710947974, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -1.379259372453352, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 1.244551457097626, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -1.577199557479029, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.09561234425855132, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.05734822654625236, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0383600177566906, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 1.379259372453352, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.244551457097626, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 1.577199557479029, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.09561234425855132, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.05734822654625236, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0383600177566906, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 8;
		at2.atomno() = 10;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 7.887211580240174, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -5.640311005121711, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 3.512323829977876, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 16.17316389050687, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.4454906298437135, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.06223766152228326, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.1688788008557868, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 5.640311005121711, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -3.512323829977876, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -16.17316389050687, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.4454906298437135, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.06223766152228326, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.1688788008557868, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 10;
		at2.atomno() = 17;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 8.045684043636173, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 20.37027447715524, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -14.94399284985846, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -17.4658018995711, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.132773070347892, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.5132517403909531, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.5939981619163894, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -20.37027447715524, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 14.94399284985846, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 17.4658018995711, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.132773070347892, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.5132517403909531, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5939981619163894, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = 2;
		at1.atomno() = 10;
		at2.atomno() = 18;
		AtomPairConstraint atom_pair_cst( at1, at2, func );
		atom_pair_cst.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ], 8.205549671860172, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -11.05052683809827, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 9.147387677649728, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -9.384204771340356, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.6086060255171939, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.4122250317625372, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.3148519362947217, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		atom_pair_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 11.05052683809827, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -9.147387677649728, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 9.384204771340356, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.6086060255171939, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.4122250317625372, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.3148519362947217, 1e-12 );
		}
		}
		std::cout.precision( before_precision );
	}


};
