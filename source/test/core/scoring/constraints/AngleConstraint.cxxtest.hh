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

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/FourPointsFunc.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.hh>

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
#include <core/scoring/constraints/FourPointsFunc.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
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

		{ // Atom1 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 1" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), 0.2175525558834286, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2460152772267948, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.0238203308330751, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01672956054012698, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01582017172013084, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0105978913123148, 1e-14 );
		}
		{ // Atom2 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 2" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), -0.2398662196545918, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1535413062411104, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.843906020387814, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.03820460893945504, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004834680065124365, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009979395600664041, 1e-14 );
		}
		{ // Atom3 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 3" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), 0.0223136637711632, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.09247397098568429, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.8677263512208889, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.02147504839932805, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01098549165500647, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0006184957116507617, 1e-14 );
		}


		//// Make sure weights are correctly applied
		weights[ angle_constraint ] = 0.5;
		{ // Atom1 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 1" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), 0.1087762779417143, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1230076386133974, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.01191016541653755, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.008364780270063488, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007910085860065419, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005298945656157401, 1e-14 );
		}
		{ // Atom2 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 2" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), -0.1199331098272959, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.07677065312055521, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.421953010193907, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01910230446972752, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002417340032562183, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004989697800332021, 1e-14 );
		}
		{ // Atom3 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		/*std::cout << "Atom 3" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
		std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;*/
		TS_ASSERT_DELTA( f1.x(), 0.0111568318855816, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.04623698549284214, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4338631756104445, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01073752419966402, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005492745827503234, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0003092478558253808, 1e-14 );
		}

		/*for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang )
			{
				at1.rsd() = at2.rsd() = at3.rsd() = ii;
				at1.atomno() = ( rsd_type.bondangle( bondang ) ).key1();
				at2.atomno() = ( rsd_type.bondangle( bondang ) ).key2();
				at3.atomno() = ( rsd_type.bondangle( bondang ) ).key3();

				std::cout <<"{\n"
				"at1.rsd() = at2.rsd() = at3.rsd() = " << ii << ";\n"
				"at1.atomno() = ( ubqstump->residue_type( " << ii << " ).bondangle( " << bondang << " )).key1();\n"
				"at2.atomno() = ( ubqstump->residue_type( " << ii << " ).bondangle( " << bondang << " )).key2();\n"
				"at3.atomno() = ( ubqstump->residue_type( " << ii << " ).bondangle( " << bondang << " )).key3();\n"
				"AngleConstraint ang_cst2( at1, at2, at3, func );\n"
				"ang_cst2.score( cfunc, weights, emap );\n";

				AngleConstraint ang_cst2( at1, at2, at3, func );
				ang_cst2.score( cfunc, weights, emap );
				std::cout << "TS_ASSERT_DELTA( emap[ angle_constraint ], " << emap[ angle_constraint ] << ", 1e-14 );\n";
				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 1\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;
				std::cout << "}\n";
				}

				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 2\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;
				std::cout << "}\n";
				}

				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 3\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-14 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-14 );" << std::endl;
				std::cout << "}\n";
				}

				std::cout  << "}\n";

			}
		}*/

		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 1 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 1 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 1 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.0787470701274538, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2313083619718075, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2310434679328141, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.07179373151825814, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006934214089804028, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.009442355675871771, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008046007442686022, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4211996216342746, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.5248243011770679, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.49956389528489, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.02526233846689022, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004301208246570523, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01678085400236346, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1898912596624669, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2937808332442537, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5713576268031483, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01832812437708619, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005141147429301253, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008734846559677427, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 2 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 2 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 2 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.1226119785767397, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3571078027582377, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3604624875263657, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.07718520469399683, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.00244597202643612, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0005720596050474271, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01398819538664041, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5592016662820729, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5418531109316609, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.7846603122002285, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01746696774303061, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01309225640958269, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0214890916098004, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2020938635238348, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1813906234052948, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.7074751075062313, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01502099571659448, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01252019680453525, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.007500896223159969, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 3 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 3 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 3 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.1664768870260259, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1257994407864308, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1294190195935522, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.005391473175738907, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009380186116240164, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008870296070824351, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005942187943954419, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2487909143855457, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3508547190226142, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.6769078754482198, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003073209381902674, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.02750498831924294, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01312685034187666, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1229914735991145, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2214356994290616, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.6822993486239582, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006306976734337512, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01863469224841856, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.007184662397922221, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 4 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 4 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 4 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.1952452586182816, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3317633141709167, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4203200520720319, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2759468273910058, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.003557102789848976, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00704790491584555, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01501192099506943, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1821120460709263, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2796458719201222, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.7088008077388908, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003023120513483506, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02322406633483865, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009939408738034979, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1496512680999906, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.14067418015191, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4328539803478844, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0005339822763654722, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01617616141899309, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005072512257034459, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 5 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 5 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 5 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.2240136302105373, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1304161094689376, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1088157223104006, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3298802193433022, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01587673854290142, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001596597332281314, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005750110096464008, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3630851059548278, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3227133305369297, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.7815073189700827, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01688271073566749, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01349900130569591, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01341787623994018, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2326689964858899, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.213897608226529, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4516270996267803, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.001005972192766063, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01509559863797721, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.007667766143476169, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 6 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 6 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 6 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.252782001802793, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0646347110455264, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.05165955362143189, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3173023751674466, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01646117926664637, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.00290063397558809, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002880905346163451, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2984457085641926, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3970327427888496, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5891370444478108, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02214267724937184, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001762671515934779, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01240496651276488, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3630804196097188, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4486922964102817, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2718346692803643, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005681497982725475, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004663305491522867, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01528587185892833, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 7 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 7 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 7 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.2876641634809608, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1087762779417143, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1230076386133974, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.01191016541653755, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.008364780270063488, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007910085860065419, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005298945656157401, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1199331098272959, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.07677065312055521, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.421953010193907, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01910230446972752, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002417340032562183, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004989697800332021, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0111568318855816, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.04623698549284214, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4338631756104445, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01073752419966402, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005492745827503234, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0003092478558253808, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 8 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 8 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 8 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.3271923948938076, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1947148004274582, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2091225374111478, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.08211134458009471, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006951218541672617, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.009214689595945747, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.006984347188852069, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2599148822670422, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2097735079113009, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.526372322686925, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0009119478171015655, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02092239597056754, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.007887830359213364, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06520008183958394, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.0006509705001530333, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4442609781068305, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.006039270724571058, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0117077063746218, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0009034831703612938, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 9 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 9 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 9 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.3640066591305239, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3168036654053773, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3468493807685167, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.07187522571084622, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.001773456277222105, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001044246117714345, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01285607883969186, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5153648746277935, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5875271681076059, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4905890662016062, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.009470961234631703, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009514335199341764, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.02134358852808888, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1985612092224158, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2406777873390887, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5624642919124521, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0112444175118538, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0105585813170561, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.008487509688397012, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 10 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 10 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 10 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.4016951433501791, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2530508119719195, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2816485803399988, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.196114423414602, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.00827525544928384, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0009040280347367525, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009379431997494801, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4395194543570223, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4965858618402317, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.378366825758009, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.006333206911169622, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.008277660798159421, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01822077544348921, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1864686423851027, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2149372815002327, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1822524023434068, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.001942048538114223, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009181688832896173, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008841343445994406, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 11 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 11 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 11 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.4371547044837797, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2762139245664825, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2429534776572204, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2723612745814353, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003171507994955575, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.006980844688582238, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009443468672664086, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2242384556967263, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1706435866475251, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5465453236652549, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.00231682278901284, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.02304967414918291, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.006246066256040151, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.05197546886975628, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.07230989100969532, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2741840490838194, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005488330783968416, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01606882946060066, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003197402416623935, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 12 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 12 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 12 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.4718550194235597, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.241287446755204, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2058280373318188, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2571139474105323, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.007714486720990676, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.002595472992684147, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009317386085405933, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1018390974982273, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.04532250377565012, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5359341824200711, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.02332854365904501, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.002166835360817779, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.004249684948399301, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1394483492569767, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1605055335561686, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2788202350095386, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01561405693805433, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004762308353501925, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005067701137006629, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 13 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 13 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 13 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.51486228934235, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1706366445539818, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.127177662054306, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4398227234584118, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009012374085023832, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.008025251687950246, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005817057382317209, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1606800316453356, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.04528226137909352, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.8830664205021932, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008306662935622846, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02458979773938751, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.002772380938757396, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.009956612908646129, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.08189540067521252, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4432436970437821, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.000705711149400982, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01656454605143729, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.003044676443559813, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 14 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 14 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 14 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.5546339815271565, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.277635354044644, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2729765691107609, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.124988268305719, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.006903101470611406, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001920362802365363, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01113969327230026, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2044963899202032, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1627954059447619, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2973317741218052, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01900356342939919, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00698777014979479, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009244162507910141, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07313896412444074, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1101811631659989, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1723435058160863, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01210046195878779, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005067407347429428, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.001895530764390107, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 15 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 15 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 15 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.589457045543112, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1435143124328669, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1015260924213561, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4185307331659952, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.007338531388797659, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0088341062583751, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004659339017296848, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4475603360340688, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3930279756261668, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4002738078674991, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.002273375229928978, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01819191556180838, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01532066061692243, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.304046023601202, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2915018832048107, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.01825692529849604, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009611906618726636, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009357809303433284, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01066132159962558, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 16 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 16 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 16 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.6242919932227435, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1271298336767355, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1646159053576974, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2970414466931842, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0006130043855969089, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01071587475222759, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006200934550134432, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4009197027818261, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.48757101785403, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5653251159087558, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01340525595521503, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009230740193652668, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01746795312493829, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2737898691050906, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3229551124963327, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2682836692155717, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01279225156961813, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001485134558574923, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01126701857480386, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 17 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 17 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 17 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.6608730741762759, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.002858418480994687, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.09447120763024303, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4433355547213056, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.009928563303592704, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007751014514810517, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.001715692960580526, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4282004073388978, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.5164866726156796, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5747382257012079, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01187611767121846, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01112059644955769, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01884161831926583, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4253419888579029, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.4220154649854364, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1314026709799022, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.001947554367625753, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.00336958193474717, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01712592535868529, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 18 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 18 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 18 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.6974065695133522, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0671091756809548, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.0109655065304372, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2783047212100854, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.001650683520853895, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01259577618185277, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -9.824862282070324e-05, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4771084399195321, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4121280524684048, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4733902981215067, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.000234423716158481, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01861890076710273, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01644566201850379, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.409999264238577, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.4011625459379674, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1950855769114209, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.001885107237012377, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.006023124585249944, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01634741339568307, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 19 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 19 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 19 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.7339322803906428, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1384008413070995, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1482234852462047, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1239293534820804, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008183377447204623, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0137179902862897, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.007268189352016743, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2929558139192929, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2404379464108583, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3767087674538579, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.02395124747755908, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.00872882565258888, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01305495573655907, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1545549726121933, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.09221446116465358, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5006381209359382, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01576787003035445, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004989164633700819, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.005786766384542326, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 20 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 20 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 20 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.7688623614038659, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06357211619556147, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.07973423762571234, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1061564401821822, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01059696986660615, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.006271749366795901, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.001635309596307204, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1471168871101961, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1147154708110828, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1710221726804217, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01467128135851903, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007640353484493411, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.007495674251793457, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2106890033057575, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.194449708436795, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.06486573249823943, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004074311491912887, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00136860411769751, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009130983848100656, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 21 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 21 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 21 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.8045049803092452, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02364680434094554, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.08571764023915533, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4217145350891641, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01083155969183359, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.006061856083001932, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.00183949024066799, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4260154896970209, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5447425440223528, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.533994728108126, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01285888549270211, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008453549274624368, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01888238180502584, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4023686853560751, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4590249037831972, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1122801930189616, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.002027325800868516, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002391693191622432, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01704289156435784, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 22 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 22 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 22 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.8401204386704891, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04202131667688577, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.003364012209463034, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.3121227536753619, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0001136315594914454, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01254180897083104, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0001504720360802536, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4916328213501309, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3936619487227824, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5216071044880206, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.002473667371746303, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01814821652768538, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01602815276335455, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.449611504673245, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3970259609322454, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2094843508126585, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.002587298931237746, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005606407556854333, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0161786247994348, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 23 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 23 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 23 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.8782722459449236, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1314329013134004, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.04889870528175048, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3342861193401439, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.006631854646195227, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.007329440613119619, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003679614505344216, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4244786989750707, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3377754249052096, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4605982454171129, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.001941809924653575, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02007145663803104, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01292967983086965, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2930457976616703, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2888767196234591, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1263121260769688, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008573664570848807, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01274201602491141, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009250065325525424, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 24 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 24 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 24 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.9164067269824531, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07761224922822495, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1440614299500861, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2700383936556047, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.002590889859350113, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008690483222542218, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005380894952799123, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3590240595727919, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4854357723150837, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5776385503319023, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.017008303661792, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.006106264398680756, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0157028809732751, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.281411810344567, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3413743423649978, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.3076001566762976, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01441741380244189, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002584218823861463, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01032198602047598, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 25 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 25 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 25 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.9529253439971129, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.09149960898995947, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1488775085912908, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2293236677507378, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.006489648576268972, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01472880522773511, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.006972623171886191, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.238866390067042, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.188084463002535, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2688351844812422, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0230374828458795, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01198582570375585, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01208373367780067, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1473667810770826, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.03920695441124421, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.49815885223198, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01654783426961053, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002742979523979258, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005111110505914478, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 26 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 26 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 26 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 0.98324537564833, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2056218966961919, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1940056218610587, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.07599380355914223, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.003995648622332789, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0009205230620459472, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008461297751276582, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1089951431749726, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.07155157543167888, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1514142243241023, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01205263984048189, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00573986517400766, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.005963659049305639, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09662675352121929, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1224540464293798, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.07542042076496012, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.008056991218149102, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004819342111961711, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.002497638701970941, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 27 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 27 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 27 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.019790728160664, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.09693429172046171, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1178465751019844, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.08534019173404372, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008845891880341702, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005290369069246265, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002742182592580317, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2133096844140565, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1421185188576697, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2370356034892502, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01847430622679316, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008428525260938232, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01157167473742519, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3102439761345182, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.259965093959654, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1516954117552065, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.009628414346451454, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.003138156191691966, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01431385733000551, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 28 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 28 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 28 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.056309418085746, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0302910927015745, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1141004155043366, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3402550824389776, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.008843487868388494, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005364561496040854, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.002586229132025531, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4101962523939961, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5524490242466188, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4228175348103526, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008928372329011627, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008535994308967875, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01981489864194176, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3799051596924218, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4383486087422825, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.08256245237137495, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -8.488446062313143e-05, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.003171432812927022, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01722866950991624, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 29 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 29 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 29 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.092814866734421, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.06664659330931519, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.003710197993768156, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2550745840378084, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 4.399974928429313e-07, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01065875646209266, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0001551523501448354, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5089998010547792, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3561625193050021, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5119837629101064, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.004765116630530874, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01626773042268801, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01605402331695698, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4423532077454639, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.352452321311234, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2569091788722982, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004764676633038031, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005608973960595363, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01589887096681215, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 30 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 30 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 30 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.129319186046143, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2169037835797199, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.02313015233485863, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5584153372040214, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0120773246062896, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01159695234714051, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005171521060817533, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5641348894927701, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3591176140003645, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.6795707072500011, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004785468376667303, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0231028769871155, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01618124444842186, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3472311059130503, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3359874616655057, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1211553700459801, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.007291856229622282, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01150592463997499, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01100972338760433, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 31 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 31 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 31 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.165792365421548, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.09341548312618954, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2369151702362725, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4067319912348417, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.002447151248760886, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01473736028972513, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009146332519271713, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4236673713681547, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.636078576465707, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.6428231430598742, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0145772583360447, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01199338074725645, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02147502220906641, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3302518882419652, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3991634062294347, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2360911518250329, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01213010708728382, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002743979542468684, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0123286896897947, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 32 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 32 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).bondangle( 32 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.202276521076142, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03289747742541677, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1026610106529422, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2038942995090342, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.007211615556553905, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01468859441977734, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.006232160444371278, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1451131984315144, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.05596728080688413, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2894046061919847, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02411614190404568, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0118213662713072, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.00980620452880567, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1122157210060975, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.04669372984605806, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4932989057010184, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01690452634749176, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002867228148470143, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003574044084434387, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 1 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 1 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 1 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.24552219270234, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2212719925717752, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1401754655867694, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4805338157410289, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.009717983249869966, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007393229938419094, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006631514498873838, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1939774050918646, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.06383025702373178, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -1.023888497711854, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01054436186257722, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.02776045448217521, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003728262311448962, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02729458747991077, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.07634520856303763, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5433546819708247, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0008263786127072506, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02036722454375612, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002903252187424877, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 2 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 2 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 2 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.28165925967871, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2729387672709551, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2803513847439175, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.183399483811536, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008114200408656293, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.001592248397070412, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009641743679989189, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2690695932804215, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2633118420453003, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1069918903071684, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01476645972705408, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01197562523022973, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.007663023539264192, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.003869173990533452, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.017039542698617, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.07640759350436753, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006652259318397793, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01038337683315931, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.001978720140724993, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 3 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 3 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 3 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.317495582945819, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3152316735163121, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2709202820203382, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2388489534522063, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.002488911483551044, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.006445106571254752, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01059537329157955, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5463272934245668, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.4820892769817167, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1739802260552124, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.007373117929833836, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001489572789068728, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01902531436395564, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2310956199082549, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2111689949613787, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.06486872739699398, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004884206446282796, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.007934679360323485, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.008429941072376098, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 4 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 4 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 4 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.354064358228394, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06456885357986743, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.006993657477817334, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4623102795785207, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01138912413754224, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005545255450122382, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.001674555683232578, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04880766180841405, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1668260795177393, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.9058424664259616, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0109376440820281, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02191512467491338, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.004625368450130125, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1133765153882813, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1738197369955567, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4435321868474409, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0004514800555141363, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.016369869224791, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006299924133362701, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 5 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 5 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 5 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.393898884991677, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1894755699137143, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.09385260845718851, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4228838945693594, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004042490975377863, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01167160559642264, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004401595651766443, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0474811903503201, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1556947708176849, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.8321056908655566, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01167335899800324, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01837367372009328, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004103985752084485, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2369567602640346, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2495473792748735, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4092217962961969, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.007630868022625373, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.006702068123670635, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008505581403850927, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 6 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 6 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 6 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.427596090627135, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1852709553711306, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1064341720802278, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3310367357407929, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01032318256648458, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.001113129353387715, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.006135454704639102, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6435878993897262, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.5751664505868617, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1509567819735335, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006087623201272284, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0009576212644539777, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02230525161360972, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4583169440185957, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.4687322785066339, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1800799537672591, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004235559365212296, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.002070750617841693, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01616979690897062, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 7 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 7 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 7 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.464635659521812, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01060157659394144, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.02275199641622491, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3336795645584113, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01263980332172756, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001772536789175597, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0005224491165747826, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3529485698629887, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3402543753726076, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1020659637935341, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01740666248919356, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01501276678977103, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01014537076311847, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3635501464569302, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3175023789563827, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2316136007648769, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.004766859167466011, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01324023000059544, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01066781987969326, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 8 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 8 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 8 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.510545579506214, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.158786070561594, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1164092698944411, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2270632562781486, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01191100073180568, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004417351804826224, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006064743044581738, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6015225140590584, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.514854751630497, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1227780791031334, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0126505110778229, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009474610017487362, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.02224761341982114, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4427364434974645, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3984454817360559, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1042851771750157, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0007395103460172072, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005057258212661135, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0161828703752394, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 9 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 9 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 9 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.554410487955501, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1552093024853834, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.113787069100248, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2219484965044249, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01164269705112032, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004317847844176458, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005928130435959591, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1201321872825993, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2228462980677579, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.6504176603806756, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02647436012991764, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004752367681984863, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003261558493177249, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2753414897679827, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3366333671680056, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4284691638762505, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01483166307879731, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0004345198378084051, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009189688929136835, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 10 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 10 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 10 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.596276974931831, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4227905430582648, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3804949515145657, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.09958698304221465, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0007061943632196352, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004829421606285998, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01545380926887297, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6917868419460096, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.7093706684528944, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5181821586815977, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01378367593532634, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004404915193314037, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02443172361062555, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2689962988877448, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3288757169383287, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4185951756393831, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01448987029854597, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0004245064129719599, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008977914341752576, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 11 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 11 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 11 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.638717943973247, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2557045109772587, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2609996438638212, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1818380053444227, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.002111004611905139, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009054058479944943, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01002712625106247, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6752374809006902, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.5978845751449675, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.008984304801265051, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004057414101568619, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004239651948259669, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02280594632514212, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4195329699234319, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3368849312811464, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1908223101456876, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.001946409489663481, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004814406531685272, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01277882007407965, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 12 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 12 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 12 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.673456043107447, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1702827082612341, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2112846864496893, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4001031730739606, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009383741822514028, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004758942333064469, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006506773215103523, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1056959574699932, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.003194443058620737, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.9395056666829725, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01440100104554534, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01890338858945392, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.001684410694494818, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2759786657312274, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2144791295083101, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5394024936090112, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005017259223031301, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01414444625638944, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.008191183909598341, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 13 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 13 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 13 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.708188171341911, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07144976900456287, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.03541780821789976, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.228467441118739, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0113897003456807, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.003801603449283446, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.002972620490330097, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2509497336617322, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2198747574919567, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.02748882654443205, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01466361856861974, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.017565091802311, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.006631389875813281, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3223995026662952, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2552925657098566, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2559562676631713, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003273918222939034, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01376348835302755, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009604010366143381, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 14 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 14 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 14 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.743617017043396, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2527122363436787, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2433155082621767, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2587916357429909, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009722490007756329, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001782383529562171, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.007818299971197656, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2194992924356891, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1933517913834451, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.0133235036334497, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0154443211538408, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01791564350351409, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.005554409846689021, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0332129439079897, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.04996371687873164, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2721151393764406, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005721831146084471, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01613325997395192, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002263890124508636, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 15 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 15 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 15 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.779001220336937, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1417554476813114, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.07343888803100397, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4382419427230745, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01155344968861593, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.002744469118592227, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004197031397198708, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06826612671031214, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.03196391391236277, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.9599333707320781, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01474766587004712, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01927657726711919, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0004069148782450215, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07348932097099929, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1054028019433667, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5216914280090026, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003194216181431169, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01653210814852695, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.003790116518953689, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 16 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 16 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 16 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.815475707038309, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3575305056442519, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3059146203182654, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.307276463539371, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01244534686047874, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002394919489537653, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01209643663167873, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.8151355483554513, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.7231085947236493, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.003187366451737021, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.004954939313511193, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005707136381082274, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.02758772653384952, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4576050427111997, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.4171939744053841, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3104638299911081, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.007490407546967544, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.003312216891544623, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0154912899021708, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 17 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 17 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 17 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.853253107684609, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3786766326811263, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3201183442008898, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1625592359711269, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.001385002996694217, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004518915408924911, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01212515534755869, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6276881269708521, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5388465331187021, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.0366501920012112, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.003058519372702061, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004952542375287198, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.02043274405402289, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2490114942897259, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2187281889178123, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1992094279723381, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.001673516376007845, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.009471457784212104, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.008307588706464203, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 18 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 18 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 18 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.899168001105266, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3176961744290899, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3028238176138757, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2191258596937211, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.009643207797540359, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00267824669833261, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01027981517083832, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4575249053206464, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3942552717372512, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.04047882000402193, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0184468927568352, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02014161367841312, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01232682889402574, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1398287308915565, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.09143145412337567, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1786470396896992, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.008803684959294841, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.0174633669800805, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002047013723187423, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 19 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 19 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 19 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.92640020152719, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.08115294986430764, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.03603137770327662, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3728701665143382, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01057087302694024, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.002026689353912541, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002496531023324655, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01746087072446213, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.09724047483127574, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.7621486923474433, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0102146270225408, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01636478194681837, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.002321955632446773, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09861382058876966, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.1332718525345523, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3892785258331053, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.0003562460043994385, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01433809259290583, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.004818486655771427, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 20 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 20 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 20 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.962381400304628, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2163018097664606, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2074767811487881, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4237879609526097, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.006955971956726697, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.007365325248177216, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.007156228059696054, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1339329592841482, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.01339013007884149, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.9322325514841099, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01111817047212132, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.02163152837477013, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.001908041560088655, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3502347690506086, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2208669112276297, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5084445905315, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.004162198515394621, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01426620312659291, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009064269619784706, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 21 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 21 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 21 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 1.998357196249527, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.002675690115215643, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.02323330340192982, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3648431104844958, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.0123260686715261, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.001363336728057076, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0001772146823876055, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4230205686704763, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3803033163559154, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1834583428465727, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01682585691076654, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01292054680680263, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01201339075020878, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.425696258785692, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3570700129539857, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.181384767637923, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.00449978823924045, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01155721007874556, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01219060543259638, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 22 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 22 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 22 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.034793981872019, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3031694758198432, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.1877358074891367, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.4659515604265967, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01470749629753806, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.0001373574355364336, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009514029806045019, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.7882088735155592, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.6509987385054734, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1913243569819886, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01114120667543544, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.005957002782579343, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02562975641924086, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4850393976957161, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.4632629310163369, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.2746272034446083, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.003566289622102619, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005819645347042909, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01611572661319584, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 23 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 23 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 23 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.079368242825755, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0214597518041886, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.02123202964765318, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.3656986375080859, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01357578623117487, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.002327323671622515, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0009317694224545026, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4203912488124043, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.282042296094405, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.9148628833757918, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02318880791552393, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.005965389383651971, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01249461995305922, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3989314970082159, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3032743257420583, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5491642458677066, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.009613021684349077, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.00829271305527449, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01156285053060473, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 24 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 24 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 24 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.118925840204084, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01522941068715523, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.02497066161519541, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.3511650071590086, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01283227148239417, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.001995455286980049, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0006984060660771632, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3278874258697021, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2880349419968531, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.2626436968967792, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01800402337832213, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01208106146300852, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.009227425107446497, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3431168365568574, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2630642803816577, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.08852131026222954, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.005171751895927973, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01008560617602848, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.009925831173523662, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 25 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 25 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 25 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.16657235686875, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4038439605134482, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3057489466493581, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.5782682802419319, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01000621402252147, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.008870972294392973, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01167838487101729, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.772970577252368, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.5881818062410258, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.6838827786174654, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.004385513150914543, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.02020016797571067, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.02233020686528791, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3691266167389191, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.282432859591667, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.1056144983755331, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005620700871606923, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01132919568131768, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0106518219942706, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 26 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 26 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 26 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.210454513233034, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3644078728273155, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.2922863131762738, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.107313811412619, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.005446223801039541, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.01062696000978378, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01045036157499866, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3869510262443044, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3272801165241034, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.4549379409338655, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02584076259576882, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01529253349164521, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01097768951084207, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02254315341698912, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.0349938033478296, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.5622517523464849, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.02039453879472929, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.004665573481861431, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.0005273279358434087, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 27 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 27 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 27 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.254343545580271, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3641613120613822, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2920594553734815, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.1076585546277012, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.005444585707788573, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.01063783896232034, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01044203112882854, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.8124097870187511, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.5921638310589493, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.7660688226582135, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.007773355711675441, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.02001500548266728, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.02371498753909065, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4482484749573697, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.3001043756854684, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -0.6584102680305128, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.01321794141946402, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.009377166520346963, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.01327295641026213, 1e-14 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 28 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 28 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).bondangle( 28 )).key3();
		AngleConstraint ang_cst2( at1, at2, at3, func );
		ang_cst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ angle_constraint ], 2.298167083397485, 1e-14 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02214017952842918, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.03466174115659808, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.561449172177452, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.02037818896933914, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004675687700183229, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.0005149340310385337, 1e-14 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4261330837584403, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), -0.2655236933160862, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), -1.218951735170436, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), -0.03358353538986884, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), 0.004681240562313574, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), -0.012760173626093, 1e-14 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		ang_cst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4482732632868695, 1e-14 );
		TS_ASSERT_DELTA( f1.y(), 0.3001854344726841, 1e-14 );
		TS_ASSERT_DELTA( f1.z(), 0.6575025629929834, 1e-14 );
		TS_ASSERT_DELTA( f2.x(), 0.01320534642052968, 1e-14 );
		TS_ASSERT_DELTA( f2.y(), -0.009356928262496799, 1e-14 );
		TS_ASSERT_DELTA( f2.z(), 0.01327510765713153, 1e-14 );
		}
		}
		std::cout.precision( before_precision );

	}


};
