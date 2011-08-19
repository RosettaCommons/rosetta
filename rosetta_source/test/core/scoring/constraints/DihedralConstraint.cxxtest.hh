// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/DihedralConstraint.cxxtest.hh
/// @brief  test suite for dihedral constraints
/// @author Ian Davis

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/conformation/Residue.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/FourPointsFunc.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/types.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <utility/stream_util.hh>
#include <basic/datacache/CacheableData.hh>
#include <utility/keys/Key2Tuple.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class DihedralConstraintTests : public CxxTest::TestSuite
{

public:
	DihedralConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_dihe_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		CircularHarmonicFuncOP func = new CircularHarmonicFunc( numeric::conversions::radians( 80 ), 10 );

		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 ), at4( 4, 1 );

		DihedralConstraint dcst( at1, at2, at3, at4, func );
		EnergyMap weights, emap;
		weights[ dihedral_constraint ] = 1.0;
		dcst.score( fourpts, weights, emap );
		//Size before_precision = std::cout.precision();
		//std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ dihedral_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ dihedral_constraint ],  0.02467401100272339, 1e-12 );

		//std::cout.precision( before_precision );
	}

	void test_dihedral_derivatives() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->total_residue() == 2 );
		AtomID at1( 1, 1), at2( 2, 1 ), at3( 3, 1 ), at4( 4, 1 );
		CircularHarmonicFuncOP func = new CircularHarmonicFunc( numeric::conversions::radians( 80 ), 10 );
		DihedralConstraint dcst( at1, at2, at3, at4, func );
		EnergyMap weights, emap;
		weights[ dihedral_constraint ] = 1.0;
		ConformationXYZ cfunc( ubqstump->conformation() );
		dcst.score( cfunc, weights, emap );
		Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		//std::cout << "Dihedral constraint func: " << emap[ dihedral_constraint ] << std::endl;
		TS_ASSERT_DELTA( emap[ dihedral_constraint ],  0.002860920447520467, 1e-12 );


		{ // Atom1 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 1" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		///Atom 1
		//f1: 0.1740739808994528 -0.1893926351891285 -0.0506199541394906
		//f2: 0.001384097747742925 0.00308827586382042 -0.006794974530963005
		TS_ASSERT_DELTA( f1.x(),  0.1740739808994528,   1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1893926351891285,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0506199541394906,   1e-12 );
		TS_ASSERT_DELTA( f2.x(),  0.001384097747742925, 1e-12 );
		TS_ASSERT_DELTA( f2.y(),  0.00308827586382042,  1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006794974530963005, 1e-12 );

		}


		{ // Atom2 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 2" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		//Atom 2
		//f1: -0.3089425719418946 0.3180238257099374 0.01152220684699815
		//f2: -0.004090960041207454 -0.004396778130442517 0.01166516855527395
		TS_ASSERT_DELTA( f1.x(), -0.3089425719418946,   1e-12 );
		TS_ASSERT_DELTA( f1.y(),  0.3180238257099374,   1e-12 );
		TS_ASSERT_DELTA( f1.z(),  0.01152220684699815,  1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004090960041207454, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004396778130442517, 1e-12 );
		TS_ASSERT_DELTA( f2.z(),  0.01166516855527395,  1e-12 );
		}
		{ // Atom3 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 3" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		//Atom 3
		//f1: 0.345548913740932 -0.3712128012163863 0.1668025760672934
		//f2: 0.008766131198070482 0.002479039606067929 -0.01264294548939172

		TS_ASSERT_DELTA( f1.x(),  0.345548913740932,    1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3712128012163863,   1e-12 );
		TS_ASSERT_DELTA( f1.z(),  0.1668025760672934,   1e-12 );
		TS_ASSERT_DELTA( f2.x(),  0.008766131198070482, 1e-12 );
		TS_ASSERT_DELTA( f2.y(),  0.002479039606067929, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01264294548939172,  1e-12 );
		}


		{ // Atom4 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 4" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		//Atom 4
		//f1: -0.2106803226984901 0.2425816106955775 -0.1277048287748009
		//f2: -0.006059268904605954 -0.001170537339445831 0.007772751465080771
		TS_ASSERT_DELTA( f1.x(), -0.2106803226984901,   1e-12 );
		TS_ASSERT_DELTA( f1.y(),  0.2425816106955775,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1277048287748009,   1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006059268904605954, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001170537339445831, 1e-12 );
		TS_ASSERT_DELTA( f2.z(),  0.007772751465080771, 1e-12 );
		}

		//// Make sure weights are correctly applied
		weights[ dihedral_constraint ] = 0.5;

		//Atom 1
		//f1: 0.08703699044972639 -0.09469631759456426 -0.0253099770697453
		//f2: 0.0006920488738714626 0.00154413793191021 -0.003397487265481502
		//Atom 2
		//f1: -0.1544712859709473 0.1590119128549687 0.005761103423499073
		//f2: -0.002045480020603727 -0.002198389065221258 0.005832584277636977
		//Atom 3
		//f1: 0.172774456870466 -0.1856064006081932 0.08340128803364671
		//f2: 0.004383065599035241 0.001239519803033965 -0.006321472744695862
		//Atom 4
		//f1: -0.105340161349245 0.1212908053477888 -0.06385241438740043
		//f2: -0.003029634452302977 -0.0005852686697229157 0.003886375732540386

		{ // Atom1 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 1" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		TS_ASSERT_DELTA( f1.x(), 0.5 * 0.1740739808994528,   1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5 *-0.1893926351891285,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5 *-0.0506199541394906,   1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5 * 0.001384097747742925, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.5 * 0.00308827586382042,  1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5 *-0.006794974530963005, 1e-12 );

		}


		{ // Atom2 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 2" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		TS_ASSERT_DELTA( f1.x(), 0.5 *-0.3089425719418946,   1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5 * 0.3180238257099374,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5 * 0.01152220684699815,  1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5 *-0.004090960041207454, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.5 *-0.004396778130442517, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5 * 0.01166516855527395,  1e-12 );
		}

		{ // Atom3 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 3" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		TS_ASSERT_DELTA( f1.x(), 0.5 * 0.345548913740932,    1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5 *-0.3712128012163863,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5 * 0.1668025760672934,   1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5 * 0.008766131198070482, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.5 * 0.002479039606067929, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5 *-0.01264294548939172,  1e-12 );
		}


		{ // Atom4 derivative vectors:
		Vector f1( 0.0 ), f2( 0.0 );
		dcst.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		//std::cout << "Atom 4" << std::endl;
		//std::cout << "f1: " << f1.x() << " " << f1.y() << " " << f1.z() << std::endl;
		//std::cout << "f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;
		TS_ASSERT_DELTA( f1.x(), 0.5 *-0.2106803226984901,   1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5 * 0.2425816106955775,   1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5 *-0.1277048287748009,   1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.5 *-0.006059268904605954, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.5 *-0.001170537339445831, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.5 * 0.007772751465080771, 1e-12 );
		}


		/*for ( Size ii = 1; ii <= ubqstump->total_residue(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe )
			{
				at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = ii;
				at1.atomno() = ( rsd_type.dihedral( dihe ) ).key1();
				at2.atomno() = ( rsd_type.dihedral( dihe ) ).key2();
				at3.atomno() = ( rsd_type.dihedral( dihe ) ).key3();
				at4.atomno() = ( rsd_type.dihedral( dihe ) ).key4();

				std::cout <<"{\n"
				"at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = " << ii << ";\n"
				"at1.atomno() = ( ubqstump->residue_type( " << ii << " ).dihedral( " << dihe << " )).key1();\n"
				"at2.atomno() = ( ubqstump->residue_type( " << ii << " ).dihedral( " << dihe << " )).key2();\n"
				"at3.atomno() = ( ubqstump->residue_type( " << ii << " ).dihedral( " << dihe << " )).key3();\n"
				"at4.atomno() = ( ubqstump->residue_type( " << ii << " ).dihedral( " << dihe << " )).key4();\n\n"
				"DihedralConstraint dcst2( at1, at2, at3, at4, func );\n"
				"dcst2.score( cfunc, weights, emap );\n";

				DihedralConstraint dcst2( at1, at2, at3, at4, func );
				dcst2.score( cfunc, weights, emap );
				std::cout << "TS_ASSERT_DELTA( emap[ dihedral_constraint ], " << emap[ dihedral_constraint ] << ", 1e-12 );\n";
				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 1\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
				std::cout << "}\n";
				}

				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 2\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
				std::cout << "}\n";
				}

				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 3\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
				std::cout << "}\n";
				}

				{
				Vector f1( 0.0 ), f2( 0.0 );

				std::cout << "{\n"
				"//Atom 4\n"
				"Vector f1( 0.0 ), f2( 0.0 );\n"
				"dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );\n";

				dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );

				std::cout << "TS_ASSERT_DELTA( f1.x(), " << f1.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.y(), " << f1.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f1.z(), " << f1.z() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.x(), " << f2.x() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.y(), " << f2.y() << ", 1e-12 );" << std::endl;
				std::cout << "TS_ASSERT_DELTA( f2.z(), " << f2.z() << ", 1e-12 );" << std::endl;
				std::cout << "}\n";
				}
				std::cout  << "}\n";

			}
		}*/
		std::cout.precision( before_precision );

		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 1 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 1 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 1 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 1 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.04672582889680656, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4358830462779327, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4990298170698553, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.05086477237590595, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01161869197634262, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00815266059184435, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01958074482810958, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5124163636972786, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5678209713020796, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05263850595783909, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.016353474987221, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01268752333503914, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02233244165759603, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2811606297026323, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2811575113997486, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.08441564481615997, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007314470545005953, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01029080502460962, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.009912806902872221, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3576939471219782, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3499486656319731, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.08618937839809322, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002579687534127565, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005755942281414835, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01266450373235866, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 2 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 2 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 2 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 2 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.04676024451672181, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01220925493611204, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0139780207313736, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.001424742206671538, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0003254441152842529, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0002283592179677338, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0005484643357587976, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02063699476758915, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02279609635434581, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.002795257463357893, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0005154725350248996, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0003583663704206452, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0008830849144440731, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01688790158894461, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01712283543370222, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.002967982497728546, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0003822484591047684, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0002568376453780916, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0006932607003151599, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.008460161757467529, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.008304759810729987, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.001597467241042166, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0001922200393641206, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0001268304929251808, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0003586401216298846, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 3 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 3 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 3 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 3 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.09290803399408741, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4470816191946403, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5118507373643892, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05217157443486949, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01191719583054133, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00836211623569778, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02008380729083546, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6422630159365896, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.7076263007624003, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1041164223721954, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01241296716827705, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00728355397032966, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02706926835919079, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2419929593379384, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.239907118674254, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.09128341418067361, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01357475855925664, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01660925734588681, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.007664952975285403, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04681156259598931, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.04413155527624288, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.03933856624334771, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01407052989699237, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01553069508051869, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0006794919069300652, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 4 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 4 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 4 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 4 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.1367729424433734, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.08017941513085999, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08395000856989994, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1089771268302979, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0160281968569982, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01799139309215634, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002066912215673654, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03091329176377576, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02332214320857311, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1053593872364921, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02067359996235995, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02232682641978587, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.00112357889879429, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.406960070489062, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4105765309932999, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.08980711799189905, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002065715571234179, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001420508953785301, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01585499484682659, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3576939471219778, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3499486656319728, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.08618937839809318, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.00257968753412756, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005755942281414828, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01266450373235865, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 5 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 5 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 5 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 5 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.1782149199624862, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07793355829306785, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08159853593229464, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.10592463480282, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01557924078205029, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01748744712134223, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.00200901719454996, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2481091098566242, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2618009507430571, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1482424777456847, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01810250123295475, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02159789988539755, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.007844952908563036, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.46375192669106, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4683861814375243, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.09775161069513574, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004146973996276691, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0002906966980980642, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01828111937618757, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2935763751275036, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2881837666267618, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05543376775227122, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006670234447181153, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004401149462153388, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01244518366217449, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 6 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 6 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 6 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 6 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.17824387356692, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.002059944827928447, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.00215681775273459, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.002799806763496399, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.000411791494895595, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.000462229327651324, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 5.310247176922126e-05, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.003716879642814319, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.00373564224698373, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.003962413672811948, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0006403661854009539, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.000717138243678032, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -7.541057199508456e-05, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.002829478427885456, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.00268423884919164, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.002147965697512841, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0005810156056832688, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0006439242855684817, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 3.932812372800888e-05, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.00117254361299954, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.001105414354942469, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0009853587881971586, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.000352440915177911, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0003890153695417753, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.702002350214114e-05, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 7 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 7 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 7 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 7 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.17824387356692, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -7.760343015677828e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 8.480946295705054e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 3.104137206271131e-28, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -6.928877692569489e-30, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -1.524353092365288e-29, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 2.910128630879185e-29, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 1.11277775742666e-27, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -1.263827291124675e-27, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -5.321378067893368e-28, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 5.543102154055591e-30, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 2.355818415473626e-29, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -4.365192946318778e-29, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -8.758101403407834e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 8.661097115711861e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 4.434481723244473e-28, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -6.928877692569489e-31, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -1.524353092365288e-29, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 3.25657251550766e-29, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 5.099653981731144e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -4.933360917109476e-28, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -1.552068603135566e-28, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 1.385775538513898e-30, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 6.928877692569489e-30, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.801508200068067e-29, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 8 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 8 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 8 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 8 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.2246005441862097, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5775911110513268, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.6423943168903578, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1743123984627442, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004533016840319302, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01011431925315844, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02225401640066583, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6353948629255902, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6893209917511848, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2033602616312911, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008992771994796946, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01547379961463717, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0243531048191948, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2526929467681757, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2578666113873952, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.02958083906580911, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01151442968038705, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01001429080948249, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01106337548690532, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3104966986424393, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3047932862482223, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0586287022343558, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007054674525909409, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004654810448003756, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01316246390543428, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 9 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 9 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 9 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 9 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.2662404788162838, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5474184139079534, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6088363746033364, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1652065187036564, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004296217239958746, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.009585958816514152, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02109149210928821, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.7099479868612921, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.7740700437069106, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1910661029180931, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002384021389117381, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009118790982232283, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02808479194571169, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1180630975833483, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1233129425372165, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.01150826445251681, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01527784391747789, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01521983217571699, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0063478479764977, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04446647536999002, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.04192072656635781, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0373678486669535, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01336564806663653, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01475266434143512, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0006454518599257791, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 10 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 10 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 10 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 10 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.2691013992638043, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.08703699044972639, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.09469631759456427, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0253099770697453, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0006920488738714626, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00154413793191021, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003397487265481502, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1544712859709473, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1590119128549687, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.005761103423499073, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002045480020603727, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002198389065221258, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005832584277636977, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.172774456870466, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1856064006081932, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0834012880336467, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004383065599035241, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001239519803033965, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006321472744695862, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.105340161349245, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1212908053477888, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.06385241438740043, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003029634452302977, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0005852686697229157, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003886375732540386, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 11 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 11 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 11 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 11 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.2962223682769191, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1697322388087691, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1183774875486927, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3609438450957351, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007326878960142225, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007114172370476128, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.005778642436852964, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1074541362873473, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.04651731473139143, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5771449065639042, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01348130904048838, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008929581965962578, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003229695995753422, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3866130071929531, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4453059819573141, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4127981282290988, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0154824602322984, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 1.340876531422443e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01451482238609848, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3243349046715313, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3734458091400129, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1965970667609297, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.009328030151952247, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001802000830172226, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01196587594499894, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 12 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 12 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 12 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 12 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.3637615810328435, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2780619354053355, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3303495113481125, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6964147843488108, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02180700848192971, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005222712355239851, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01118445448517348, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1407093831877251, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2266441711671725, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7261912959403833, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02088401636057385, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007441779797804006, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006369138683865899, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6491747123901569, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.6930277904889015, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2804668015512755, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0156432472457989, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005062745209733665, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02369826791751144, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5118221601725463, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5893224503079616, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.310243313142848, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01472025512444304, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002843677767169487, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01888295211620386, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 13 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 13 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 13 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 13 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.3834356056298382, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2059324220440643, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2272514355914404, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.03000759264951343, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004695626942486976, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003098265309270058, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008761002356221704, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2954817565624215, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3153486151294141, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.08896218804195501, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004982111108669475, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.008207286135938949, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01254503064418841, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.05969496722465919, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1268297664816752, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4539502803832385, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00355266424968924, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01455718357163787, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004534329987040427, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02985436729369806, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.03873258694370157, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.33498049969177, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003839148415871739, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009448162744968982, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0007503016990737214, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 14 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 14 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 14 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 14 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.4596785473939113, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4053945863449546, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4473627844195236, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05907236698626498, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.009243720454919009, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006099185213145808, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01724673993009078, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4904095379946661, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.489225793216468, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1577768620541205, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01152238536864684, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005141305007207106, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01987255053811627, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3710426727748548, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4213706984824365, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3195613536642522, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004319063450380418, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01700460545992818, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01740725394285101, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4560576244245669, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4632337072793803, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2208568585963964, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006597728364108251, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01604672525398948, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0200330645508765, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 15 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 15 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 15 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 15 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.4647956969124502, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1050248756088816, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1158975042299257, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.01530377613226342, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002394754699844614, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001580105383501848, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004468083138819296, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1787979195018171, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1776259042395632, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.06414463377716689, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004271465182544466, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001690630888100039, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.007224754750946265, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1931082238751507, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1634356811296115, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2145941275549606, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006267154830360298, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002336226321105297, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.007418944293787679, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1193351799822151, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.101707281119974, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.165753269910057, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004390444347660448, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00244675182570349, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.00466227268166071, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 16 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 16 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 16 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 16 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.469834479607913, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07506771856050341, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.05435399108413634, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1620961600477259, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003071708247002642, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.002982533498897281, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002422628205859721, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09773866723086605, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0689815774282284, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2864795940543134, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005013203255540625, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.006056463097588502, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003168701023394305, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03777954372184288, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.004974067433114892, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2939092079141731, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003884397974369123, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007855423156732612, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0003663627260018334, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01510859505148029, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.01960165377720695, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1695257739075856, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001942902965831141, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00478149355804139, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0003797100915327524, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 17 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 17 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 17 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 17 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.4897010376737062, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1490567741045402, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1079269588322224, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3218631280498215, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006099278505706802, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005922210379233014, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004810445183958424, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2406626386346112, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2041297616183596, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.398907542001874, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009079086645015033, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006402981538571039, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008753998548065647, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3244049199200923, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3326649595487297, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.03569412018060714, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006347682830192609, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00771043600652275, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01416962624227869, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2327990553900213, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2364621567625924, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1127385341326596, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003367874690884377, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008191207165860775, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01022607287817147, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 18 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 18 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 18 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 18 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.5696671281783037, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2990496748579374, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2165317352470395, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6457476647425044, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01223686253078563, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01188161423047992, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.009651101580769279, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2782697416388529, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2083953508741215, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7083328578895298, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01312998884912219, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01426408479657305, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.009354708717860627, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.492525004145588, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4101964353388989, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5926556845360934, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01646278231225011, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007289808084075759, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01872686874786158, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4717450709265026, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4020600509659807, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6552408776831171, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01735590863058665, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009672278650168852, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01843047588495289, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 19 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 19 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 19 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 19 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.64582670445653, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2142028029188999, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2685300920280959, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7280393769912235, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004076157953574633, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02409943309253098, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01008811841374711, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2333155226014475, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3094469474624789, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6107353128833071, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001823712592289572, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02148744783276677, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01158394716649629, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07785131088812712, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.03528966967596218, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7763795836011461, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009805985758205189, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02120128456263944, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -1.960705462442316e-05, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.05873859120558015, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.0762065251103449, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.6590755194932305, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007553540396920145, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01858929930287523, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001476221698124755, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 20 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 20 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 20 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 20 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.6499274525391256, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04970435416016665, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.06231064493530462, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1689367577987276, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0009458456928500657, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005592115234570102, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002340881648670256, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07530638982018735, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.102316860270479, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2189235520888949, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0008208472965106125, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007540674626615038, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003806594542518684, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1313693418988593, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1474377738795369, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1012071527428668, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001655120729559499, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005670060389228543, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006111711907522284, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1057673062388387, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1074315585443626, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.05122035845269956, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001530122333220046, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003721500997183606, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004645999013673853, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 21 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 21 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 21 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 21 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.6696817288594681, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1090921377074644, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1367606837022501, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3707858668929453, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002075961559462412, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01227368940917925, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005137815136893346, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06360925076446462, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1245926080082538, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5262179334040558, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00432908259278255, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01584569243408485, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004275083959474828, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2799516547928777, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2120016946579553, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4811027108985517, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01087942952560903, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00837936063640487, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01002312571666545, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2344687678498779, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.199833618963959, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3256706443874414, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008626308492288893, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004807357611499261, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.009160394539246927, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 22 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 22 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 22 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 22 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.7599584508868655, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01666240851476511, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.06553921185077531, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7400432135160346, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008217794292423107, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02022403133947636, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001606039765160617, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04250745100956666, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.07390281914454412, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.796423472202933, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008549701511204882, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0232441421871677, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001700579735989101, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02314419667844755, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.07885219943019844, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4961122898760119, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0080965762323237, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01162905394787378, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001470610067978145, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04898923917325004, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08721580672396768, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.552492548562909, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008428483451105487, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01464916479556511, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001565150038806561, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 23 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 23 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 23 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 23 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.7679966068859276, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.004971965454441485, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01955651831154098, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2208245758189834, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002452141855568388, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006034732918732355, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0004792328926374638, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02370477064930656, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02295768874367198, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3196660719342132, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002524236766654155, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01022869788068883, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.000547417520792885, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07295688784310754, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03366248146359393, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1972782797626338, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004074556831833133, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01177658512215539, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0005026508377297589, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0542240826482425, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.03706365189572485, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.09843678364740402, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004146651742918902, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007582620160198912, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0005708354658851809, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 24 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 24 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 24 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 24 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.7817195780658805, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.006496418028164821, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02555273549089892, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2885315210971116, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003203992204712832, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007885040249704729, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0006261703207608712, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.008518669650851905, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.05813250231104618, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4549897204086111, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007037807945674242, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01114562992673765, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001292272264944285, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04448017869440579, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1398491643465727, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4575271241169707, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01499792409820362, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00333998860235983, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002478989548629167, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04650243031709282, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1072693975264254, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2910689248054711, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01116410835724221, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 7.939892532690523e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001812887604445755, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 25 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 25 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 25 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 25 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.7957642409611763, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.05709017401929917, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.0483441564019613, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.04414105954275849, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006971838566529678, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009153648140097871, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001008175172883313, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07474508333356902, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.05308077552763228, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1524661105905984, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007461323470569613, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01346383555823392, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001029563801510186, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03697763835242471, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.02966373959824091, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3262435978888376, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002834945209118715, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01008822861160992, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0005959504159123345, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0193227290381549, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03440035872391185, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2179185468409978, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00332443011315865, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005778041193473869, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0006173390445392055, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 26 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 26 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 26 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 26 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.8851671932300151, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1440395106222409, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1219731547363928, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1113686676228531, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01759007101503534, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02309481484610945, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002543643648283769, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1518286452262532, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1369436462292417, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1111436066698103, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0187351905002196, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0229880195251383, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002730922909124692, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1730489624733174, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1086373217188464, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3285131530299001, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01268402383654512, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02539494283147597, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001716465244429623, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1808380970773291, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1236078132116943, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3282880920768544, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01382914332172938, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02528814751050482, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001903744505270534, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 27 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 27 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 27 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 27 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.8936684785117136, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04441687521827545, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0376123632363219, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.03434230088439378, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.00542417831038362, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00712165366834027, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0007843729962529745, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06449201353520177, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.06437381668826171, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.004906244092089255, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008720666501345941, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.008835474141668374, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001296276864640428, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.05667621576224454, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1111909230432138, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2585302393132275, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0120835212945021, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001651327261787356, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001938789283999528, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03660107744531822, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08442946959127394, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2290941825209229, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008787033103539781, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -6.249321154075113e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001426885415612074, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 28 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 28 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 28 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 28 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.9019790994155303, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03018793935484534, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.06408408769718865, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2356785934951146, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00867640405617256, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001148151501119361, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001423553136404041, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.03719840135048751, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.09141833573477715, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3673230853296282, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01170217165940036, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003033332846597135, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001939993284072366, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.00785332729256074, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.05379633227533566, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2992758525931314, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005583047613094621, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.006329873740054399, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0009913211781527988, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01486378928820289, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.02646208423774717, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1676313607586178, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002557280009866822, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004444692394576624, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0004748810304844733, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 29 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 29 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 29 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 29 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 0.9162934220958688, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.03961887454601376, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08410442995219364, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3093063265798968, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01138697675821492, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00150684250912158, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001868281649099167, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06352632226480362, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1187697157601613, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3726274385856577, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01539397268553845, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0004131332363228341, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002492717005042668, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09626758793776076, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08412543187627014, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.06803933006271393, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009540557527679665, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01203871692505149, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001386195303481352, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07236014021897093, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.04946014606830249, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1313604420684748, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005533561600356135, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01011874117960707, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.000761759947537851, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 30 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 30 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 30 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 30 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.007422059812496, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.09996412568108885, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2122075879846853, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7804244027244758, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0287309821096461, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00380198063901322, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004713943637175127, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1024230266859375, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2271690759099494, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.8440448239249786, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02920446515508431, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004676558253682847, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004802563816464062, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1173748114627828, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2614647128583217, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6864451045270845, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02829569184965265, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0006699717925010412, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004583072255835863, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1198337124676273, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2764262007835832, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7500655257275789, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02876917489509043, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.000204605822168529, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004671692435124752, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 31 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 31 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 31 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 31 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.103400253858269, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01842023769572903, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08891592964372988, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7330220165939026, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01069698676406998, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01859194751130306, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001986405899236957, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.05364442050610162, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1137291211520527, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.8156121824953717, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01187196531763145, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02052913362123035, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.00208174382563535, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02547713927599356, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.09553377778273978, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4909488016211341, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00721784193498612, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01328944950692857, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002211434976957038, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06070132208636442, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1203469692910649, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.5735389675226044, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.00839282048854775, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01522663561685586, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.00230677290335547, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 32 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 32 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 32 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 32 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.114049980099187, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1801044571066397, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2190642606533986, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1125476259078291, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001217058247819599, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005406250270771218, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008575201793002809, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2249940759838366, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2630732732893191, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1627827747476167, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00102523964951328, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007425954815387403, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01058402767522267, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02466962194821916, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08409735446832303, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2412846311798217, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002987520547967985, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007091793830800532, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002777226625704975, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02021999692897769, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.04008834183240254, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1910494823400342, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002795701949661667, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005072089286184346, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0007684007434851167, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 33 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 33 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 33 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 33 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.126598859806305, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2335180387750115, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.17607734577095, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2511211984260727, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006546129899383842, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003211856811709725, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008339298518795795, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2741380290520523, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2107969166371144, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3330366386751107, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008480451290965347, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004820440168094985, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0100317735007664, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06256897020591234, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.008796667704581535, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2893012904067669, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004969079813937673, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007114379738597333, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0008583693905936678, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02194897992887151, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.04351623857074601, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2073858501577289, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003034758422356169, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005505796382212072, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0008341055913769379, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 34 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 34 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 34 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 34 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.225290519568289, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.000817202970129159, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.003985429500876811, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.02333391591229717, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0003293119940322022, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0005974527567022098, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -9.051164452006459e-05, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.001105994828825608, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.00334618409743215, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.01806193784641493, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0002680974470058229, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0004863598459468623, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 7.368728668388258e-05, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.00467427618922256, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0092639429801601, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.04415465071339458, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0006460588474320949, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001172314611879139, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0001775667592249337, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.004963068047927003, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.008624697576855323, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.03888267264746038, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0005848443004107116, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001061221701121918, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0001607424013892513, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 35 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 35 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 35 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 35 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.236240769593775, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.006808614068041733, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.03320503272598105, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1944090194847976, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002743698147579833, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004977741629628168, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0007541074601154176, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.03160587383698724, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02666209557698244, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2363017799539341, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002236152321768106, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007635846488787799, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0005624677087030321, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1300982160006803, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08462667499957122, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1510268915195426, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005375099750115148, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01145102809859305, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001786248410999053, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1053009562317348, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.07808373785057263, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1091341310504061, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.00486755392430342, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008792923239433417, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001594608659586666, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 1;
		at1.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 36 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 36 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 36 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 1 ).dihedral( 36 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.247213470659619, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.006815590264819374, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03323905504535316, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1946082135589102, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002746509377293352, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004982841890087861, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0007548801286762552, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0009088551697519837, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.06623042743625296, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2650856888097933, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.005206402068528867, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005864835004175097, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001483152094595005, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.08089586277891261, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1827128160517388, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3155542024709746, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01221465708496704, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0009243615554106879, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003666590136388509, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07317141734434124, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.149721443660839, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2450767272200917, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009754764393731526, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 4.236844132345162e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002938318170469759, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 1 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 1 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 1 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 1 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.271331272466099, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3823622441578476, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.374289630421332, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2462734800970666, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01102951615861789, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.002374161476487403, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01351605754495932, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5636696976094566, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5586975607701458, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2058900478949982, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01220798456950513, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00505508576419818, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0197047255379798, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3661583557709169, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2839757565735188, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4079249412337181, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006243277994324882, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.008444646948246389, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01148275807059207, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1848509023193079, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.09956782622470498, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3675415090316496, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007421746405212119, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005763722660535616, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005294090077571586, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 2 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 2 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 2 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 2 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.275290267313737, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1549168327103096, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.151646154784174, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.09977948424592113, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004468688359604072, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0009619087184989586, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005476128612559838, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2185390228678254, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2277260453808274, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1745632636591146, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007226401645737516, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0009916047101961689, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.007753272406265979, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.106118684758921, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1278202136717509, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2206730977913936, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006799968520156584, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0008689464196621733, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003773330963880093, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04249649460140526, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.05174032307509747, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1458893183782, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004042255234023137, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0008986424113593836, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001496187170173951, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 3 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 3 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 3 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 3 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.349867585926383, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6723718507952957, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6581764161735241, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4330640855489518, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01939505352909549, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004174887480192031, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0237675575081319, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.7268780851133564, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.7566367597863302, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5738086089030521, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02385655560467166, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003367683320246057, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0257798405781777, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6698331808586643, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5977850568705073, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1633403721437016, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003386875960443644, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009744167575111212, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02177217930639471, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.7243394151767253, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6962454004833134, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.02259584878960062, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007848378036019817, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008936963415165243, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02378446237644052, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 4 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 4 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 4 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 4 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.354443052486783, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07344326981455877, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.05104848787638519, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1585365261630411, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003287656613477219, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002553191646355754, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002345155615078564, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.136868080324214, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08565340210616007, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3050694135921692, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004916813519967651, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006047618935910166, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003903878629683202, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1476915318861162, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.07263227294701403, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3325715305091564, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002923595819941291, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.009465424597257233, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003365548597791137, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.08426672137646099, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03802735871723922, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1860386430800283, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001294438913450859, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005970997307702823, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001806825583186499, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 5 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 5 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 5 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 5 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.41521371840691, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2676586658959094, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1860425086668489, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.5777748620903564, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01198162589041941, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009304921629520011, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008546749414005518, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2638213140225306, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1772469549403009, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4976188820375743, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01390140696005391, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003507864083793297, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008619539637638408, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3404588395467144, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.15859292419203, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.8163508075316788, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002899598227689301, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02802792569260856, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.006654277171358997, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3366214876733356, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1497973704654823, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7361948274788963, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004819379297323798, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02223086814688184, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006727067394991888, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 6 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 6 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 6 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 6 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.43607382794963, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1111034213431387, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.08528599888394917, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1373237412780822, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008562783760626844, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004519771949923163, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.00412078276693862, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07102228141331111, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.06713836270180829, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.01063985839450353, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.010817208552583, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01208812170573612, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004070906578762637, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2200082880543968, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.09934403367974322, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.5451951751223048, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0005094741507614417, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0203176830570912, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003907828283984836, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1799271481245691, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0811963974976023, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3972315754497191, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002763898942717597, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01274933330127824, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003857952095808852, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 7 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 7 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 7 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 7 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.464881751999845, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1305645036984769, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1002248534032148, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1613776237422602, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01006265692332493, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005311463628554214, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004842586756648101, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2452513615511634, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1717855050435484, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4103673462749921, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01547809520159123, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00144590147375743, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008645032407424619, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3464534414534627, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1746973379675468, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7558657428625009, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008733619427477362, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01917169248086228, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008434085296141151, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2317665836007762, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1031366863272132, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5068760203297689, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003318181149211061, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0153061303260655, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004631639645364632, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 8 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 8 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 8 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 8 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.539136372085825, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1703636481347019, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.07840366419221863, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.593809443372796, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004057300279694098, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02549164271311114, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002201752990545393, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.05694659132607856, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01036870560116693, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3150604853395963, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007158837576512984, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01947177320847595, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0006531263788484317, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4528860749239646, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.221228412832492, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 1.028206804099693, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008316192617231077, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.03007406978803955, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0101336876179053, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3394690181153417, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1531934542414403, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.749457846066493, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.005214655320412191, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02405420028340436, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.007278808248511456, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 9 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 9 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 9 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 9 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.540872115664444, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02604701588060271, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01198719039344276, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.09078793610567534, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.000620323443261831, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003897434868535236, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.000336627536084201, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04841980954914289, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02438565548367761, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1519691581843713, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0004155372751293663, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006109067614372465, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0008478916642541271, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07926301468012326, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0377147500714396, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1856007495792759, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00101927830070085, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005968728024861068, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001648162290869174, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.05689022101158305, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02531628498120478, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1244195275005799, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0008144921325683871, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003757095279023839, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001136898162699247, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 10 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 10 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 10 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 10 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.624772248903357, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1789819633960427, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2436320451183561, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.6298870125543158, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01873292755973722, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004164557213609771, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006933744717248228, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1953867841995509, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2555459260023333, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.556711945293932, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01454230338006325, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005016005918104889, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.00740633249261999, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09592775975309234, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.160447848189863, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7823562567120639, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0258325429130698, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0001330472812192752, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003194714810436355, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1123325805566018, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1723617290738397, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7091811894516792, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02164191873339583, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0009844959857143937, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003667302585808097, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 11 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 11 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 11 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 11 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.640777149220286, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.07817257708767195, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1064092965943124, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2751109112400904, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008181836850794736, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001818921659142604, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003028398410268583, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07203528251411287, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1186445236560643, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3871306786861752, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01239495711938261, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00102110346832672, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0026193288940302, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1874185171316493, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1501076995126292, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1668889134504766, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01205608247476087, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007353378259028692, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006925171851337655, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1812812225580902, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1623429265743811, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05486914600439166, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007842962206172992, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006555560068212807, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.006516102335099271, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 12 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 12 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 12 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 12 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.646810098360608, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.04799464400397108, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.06533079116771753, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1689064213882557, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005023300517674665, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001116740687770276, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001859308072190373, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.09379702635599954, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.123429954180552, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2728605804554166, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007197957075564621, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.002382445741730194, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003552041620111774, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1945553122800938, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2016402178959976, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3022382533372768, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00502045192294285, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00564402307127505, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006997187186617613, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1487529299280654, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1435410548831631, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.198284094270116, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002845795365052895, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004378318017315131, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005304453638696211, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 13 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 13 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 13 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 13 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.652761816569984, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.05130334883812086, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03331392887127528, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.07006632898953065, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004527301542370062, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002389686706108076, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002178733773742695, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.06364362172809854, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.03518173904215518, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1764551038114569, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.009099236862622815, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003263085630492355, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002631305934177239, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01757862496089739, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.04403937846494031, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2952735975495189, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0103360898751741, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0006111861746616491, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0005241849914921459, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02991889785087498, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.04590718863582018, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1888848227275925, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005764154554921351, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0002622127497226281, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0009767571519266913, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 14 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 14 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 14 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 14 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.736456444410836, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1923859381901031, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1249261813477635, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2627465212986473, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0169772378299636, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008961249691236692, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008170182855909354, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2082422791668928, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1480102680196833, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3324523364371039, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02030192914446183, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009561636844108849, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008459858019545458, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3986919469828239, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3481566922907078, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.05576724940157218, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01461035186990611, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01439066467558285, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01461114594219822, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.414548287959613, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3712407789626276, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1254730645400285, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01793504318440435, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01499105182845505, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01490082110583433, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 15 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 15 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 15 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 15 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.75463712546431, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0896663741025409, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.05822503358255669, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1224597187016122, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007912674766204081, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004176619012802905, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003807922134726066, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1831089984358587, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1568015154597029, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.09387886809819661, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01017748384540277, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007503981027159673, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.007317443109463047, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3516720745215558, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3477583144625191, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3727945253808177, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002675383848908443, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01092795620034886, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01271785150867617, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.258229450188238, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.249181832585373, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3442136747774023, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004940192928107135, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007600594185992094, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.00920833053393919, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 16 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 16 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 16 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 16 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.7713451108242, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2398282300914729, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2712660233198603, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2734143511342302, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001543076939675781, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008205158697289038, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.009494213347217324, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3259153220126657, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3643428267944996, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4676159005549251, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006429806909223431, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01046616291405407, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01263610575278759, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1362157932471924, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1699936658307188, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5106754687047558, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01454449154621958, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002700338067905887, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004778435562374255, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0501287013259996, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.07691686235607953, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.316473919284061, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009657761576671928, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0004393338511408541, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001636543156803993, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 17 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 17 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 17 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 17 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.778216440847081, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1538008243642307, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1739617475085746, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1753394610026394, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0009895686812559589, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005261933389580815, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.006088598656359031, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2177241807270679, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.229044128253608, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2055401214792666, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0001635718562565872, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007478327845022125, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.008506760713516564, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1827042910709527, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1614543716751848, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.005751261940209249, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.00596494264952419, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006511795498359554, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.006687709201064322, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1187809347081155, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1063719909301516, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.03595192241683638, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.00513894582452482, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004295401042918244, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004269547143906794, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 18 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 18 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 18 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 18 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.860629083630127, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5326417196900848, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6024628589364649, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6072341447065663, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003427065923179858, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01822308340092896, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02108598359099592, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5708082210491795, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.6062531548410772, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5818040203880058, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002089308682876342, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01941043549729058, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02227594151192647, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5116244871614056, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5267375968612839, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7582883493590811, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01185581996105866, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0149949159904644, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01841530841560649, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5497909885204997, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5305278927658962, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7328582250405209, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01051806272075513, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01618226808682603, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01960526633653703, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 19 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 19 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 19 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 19 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.957130087879028, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1090705457067363, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1943740402541975, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6958401410957626, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02300875082279306, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001046673499720158, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003898917301562378, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1372135312515826, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1996016996081526, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.8087651321301205, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02586978771907299, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001078920135245103, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.004655293692598776, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.06438806872894465, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1221950686563076, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.6230230484776903, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01871701192819376, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001211821672155681, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002172039840528737, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0925310542737863, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1274227280102564, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7359480395120512, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02157804882447373, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001244068307680798, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002928416231565086, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 20 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 20 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 20 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 20 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.96494372618222, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03103617710087719, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.05530940638513698, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1980022902884829, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006547172389934115, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0002978324156500903, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.00110944239885191, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0254003733963676, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01769375086033192, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4173698392228509, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01244689847385427, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.001575608574910867, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0006907002284148704, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2567574729221457, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1624834756118509, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5211060913586599, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01297775082670642, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.004155567434691584, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.007690076191507069, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2003209224249008, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1248678200870459, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3017385424242919, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007078024742786262, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.002877791275430806, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.005889933564240291, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 21 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 21 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 21 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 21 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.978172288478057, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.04038293621389009, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.07196621616184741, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2576320476885646, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008518898579057865, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0003875267048687649, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001443558640621223, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.07991627099633616, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1088292930388146, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3940047690151908, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01272567054435849, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0003872088221999452, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.002688109082686314, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3550490928431386, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3175265018722681, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1385336797457654, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007098934689035426, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003177193436297377, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01091162240597226, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3155157580606925, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2806634249953008, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.002160958419139237, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002892162723734804, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.00317687555362856, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.009667071963907169, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 22 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 22 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 22 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 22 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 1.988060494001401, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2100610711435875, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1736480755266826, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.265475682907573, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006022177597127646, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003087921684963295, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.006784940582552897, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2511717624140653, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1837126255180101, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4048215113591185, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.009945516509423276, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003927134239816913, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.007952880371745958, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01149101304495027, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.03072413691928665, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3749266752508567, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01083058608387149, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.001237445377912597, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0002305383171129859, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0296196782255276, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.04078868691061416, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.235580846799311, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.006907247171575862, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0003982328230589778, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0009374014720800775, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 23 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 23 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 23 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 23 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.077821515526495, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6328936330466437, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.523184809033547, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7998524835002094, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01814423699521619, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009303608532130407, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02044236788815691, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5383301799032096, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4209736624725181, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5850832093394384, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01272791808509737, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.007532335873636431, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01713044792988911, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.7735224581528818, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5254326945537661, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 1.237468740883737, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02940626753166653, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01152513300088471, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02327500411313211, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6789590050094484, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4232215479927361, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -1.022699466722966, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.02398994862154776, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009753860342390737, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0199630841548643, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 24 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 24 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 24 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 24 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.087074255550435, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.2031992091721034, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1679756816871622, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2568036453759161, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.005825456942442386, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.002987051531985321, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006563303486689596, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2673155420693253, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2069225057326104, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.3125117670611147, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.006988293131141697, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003809101325831551, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008499738049548388, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3279923443317454, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2736747045117742, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0539008424164129, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001255972838937095, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003478973349956197, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01002131884174681, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2638760114345234, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.234727880466326, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.00180727926878581, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002418809027636407, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002656923556109967, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008084884278888016, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 25 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 25 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 25 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 25 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.098434323375808, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2787463288098538, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2424768886161511, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.04594036532196344, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003958468751961525, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.002837597574061951, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.009041216770993808, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3542213474863101, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3276526507747549, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1406791164219552, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007601893278398318, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003310618416758212, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01143041001013757, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1072227032054064, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1288948854565824, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.3472447489133507, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01104691832866454, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -4.617724435313622e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.00339394170352522, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03174768452895006, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.04371912329797866, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.252505997813359, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.007403493802227752, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0004268435983431243, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001004748464381455, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 26 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 26 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 26 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 26 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.11348814528689, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.320879580416185, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2791279175298043, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0528843741606142, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004556801869438978, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003266508021262937, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01040782081799999, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4972135579040479, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4314543642844179, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.09794281204897851, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001691989338686472, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005571873962207693, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01595552819569572, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.454384071653302, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3256458819546603, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.5696472930602373, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01268927530984557, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006299807087859869, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01372307199932758, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.278050094165439, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1733194352000468, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4188201068506448, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.009824462779093067, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003994441146915113, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.008175364621631848, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 27 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 27 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 27 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 27 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.204821591823671, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.7903757688749605, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.6875350003448147, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1302623489922733, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01122410399719858, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.00804591175759246, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02563606375535091, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.9278840927674736, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.8372724773261749, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1504183113516845, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01328459314735264, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.009354939694561037, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02987621046022205, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6915396346358875, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5877326289687329, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.01447783656036213, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005538946188929867, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007038517708367858, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02116101336993791, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.8290479585284004, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.7374701059500932, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.005678125799051766, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.007599435339083868, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008347545645336371, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02540116007480904, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 28 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 28 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 28 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 28 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.239966303052928, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.05232222993987351, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0871108973419952, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4225436007173263, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01333396294865145, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0007687609225071095, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.001809588709690748, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.06195365241963117, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.007972888178324594, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4147624974687875, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.01682402132125136, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004313967678424518, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002430101066232811, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.4815342615273111, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.4287435417120803, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4680046149595618, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0004600086819126789, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01803118859329509, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01604523861046898, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.3672583791678065, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.3336597561917606, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.4602235117110229, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003030049690687233, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01294845999236347, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01180554883454542, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 29 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 29 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 29 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 29 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.256546835312256, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.03593813116434007, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.05983313131182563, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2902289784415202, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.009158587276187579, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0005280323659528329, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.001242937016964079, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1198856250123947, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.156430962896432, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.5022017553381199, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01297765867210679, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003452441591691062, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.004173429227111822, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3022333218733532, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.2809647011106358, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4639431003810064, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005419586152852934, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01073211887074626, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01002994999416336, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.2182858280252985, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1843668695260294, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2519703234844067, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001600514756933722, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.007807709645008028, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.007099457784015614, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 30 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 30 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 30 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 30 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.296384062088214, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1695574368024567, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.102593400188936, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2785963375214965, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01596222312888425, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008776105376565327, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006483011099400763, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0698711920374105, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.002903946753018077, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.5241036902549023, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.02094002753057431, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.005076646760828606, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002763503482060333, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.4906948114317409, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4549266118597101, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7354930028040799, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.008203804007659932, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01748528099405128, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01628850805504076, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.3910085666666948, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.3552371584237921, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.4899856500706742, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.003225999605969871, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.01378582237831456, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01256900043770033, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 31 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 31 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 31 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 31 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.309038340168536, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0955632653037506, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.05782206021486105, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1570180359958961, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.008996374281597589, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.004946248906861472, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.003653851587653976, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1812888150805931, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1430810593128194, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1055355910572904, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.01055970758016939, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.008465923081122427, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006661673329871916, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.276423047116067, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.2463244284432223, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2716071801397494, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0001651016251606156, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01034059655261067, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.009210005567547854, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1906974973392245, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1610654293452639, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.2201247352011436, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001398231673411181, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.006820922378349716, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.006202183825329916, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 32 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 32 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 32 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 32 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.309102183476557, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.01410796851422845, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.0146702571208693, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.018188435276893, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 7.065092011840906e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0005717002519390091, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.000515917422564881, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.02338536674614972, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.02310906458371391, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.02943770093512529, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0001402059532328575, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0009305395407500062, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0008418682571662486, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.02493046416151039, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.02265985307535013, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.03086463559882687, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0001986999900946845, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0009107207691093737, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0008291198423741983, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.01565306592958907, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.01422104561250557, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.01961536994059449, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0001291449569802371, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0005518814802983772, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0005031690077728297, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 33 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 33 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 33 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 33 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.404116819677976, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5442545215949031, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5659463843387796, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7016699909842836, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.002725557736431206, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02205494340315147, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0199029640388931, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5133019362926348, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5039162741503624, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.6154180501403951, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.002005068026962191, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0205403680589927, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01849122078887166, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.5534950164044697, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5033758131520947, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.6894298319156392, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004551873907283928, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02020502179297617, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.01840674454500721, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.5225424311022029, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.4413457029636787, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.6031778910717507, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.003831384197815028, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.01869044644881744, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.01699500129498571, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 34 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 34 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 34 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 34 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.502789287708672, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.08934622946156306, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.08887572053107863, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1178844645025219, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.0007246264081085397, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.003534908110718413, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.003214250022537387, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.08294551755362209, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.07531237451756595, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.1015607514767248, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.0006217644514437089, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.003035722295557063, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.002758940088583815, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.1734404485245471, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.1464827576133825, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.2000749500417038, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.001272410659301226, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.006197681399374691, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.005640590875625571, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.1670397366165222, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.1329194115999408, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.1837512370161304, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.001169548702643389, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.005698495584211343, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.005185280941671, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 35 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 35 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 35 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 35 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.502789303117253, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0001167951911480052, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.000116180132395134, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0001541009469334516, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -9.472462392195264e-07, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 4.620903097650392e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 4.201732384704797e-06, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0002631763605362962, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.0002389370730082375, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0003219502105587788, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 1.974208132431227e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -9.619429752505381e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -8.752916469335565e-06, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.0003487702495161305, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.0002945446229956014, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.0004020288075233909, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -2.560329635449036e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 1.244900780713965e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 1.134186165691396e-05, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.0002023890801277859, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.0001717876823825587, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0002341795438979436, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 1.533367742237928e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -7.450481152285842e-06, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -6.790677572283934e-06, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 36 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 36 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 36 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 36 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.502793868581902, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.004498293445777493, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.004069818008192692, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.005336405533493951, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -3.34159647093477e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0001521913452014534, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0001442367690381247, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.008550681476986382, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.007746207605745506, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.0102173823857838, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 6.533380956864527e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0003020228058914985, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0002836519024126321, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.009337152335255226, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.007881663001152738, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.01069445571591259, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -6.891975446938616e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.0003301191375035663, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.0003034658447923106, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.005284764304046054, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.00420527340359981, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.005813478863622506, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 3.700190961008662e-05, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0001802876768135159, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.0001640507114177969, 1e-12 );
		}
		}
		{
		at1.rsd() = at2.rsd() = at3.rsd() = at4.rsd() = 2;
		at1.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 37 )).key1();
		at2.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 37 )).key2();
		at3.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 37 )).key3();
		at4.atomno() = ( ubqstump->residue_type( 2 ).dihedral( 37 )).key4();

		DihedralConstraint dcst2( at1, at2, at3, at4, func );
		dcst2.score( cfunc, weights, emap );
		TS_ASSERT_DELTA( emap[ dihedral_constraint ], 2.600206007131707, 1e-12 );
		{
		//Atom 1
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at1, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.6570702889789841, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.5944824469478042, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.7794941722384481, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.004881103879226203, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.02223074425379933, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02106881124045813, 1e-12 );
		}
		{
		//Atom 2
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at2, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.7019314902564115, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.634837808139472, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.823635158867663, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.005437438185784555, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02413877759509149, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02323955885379499, 1e-12 );
		}
		{
		//Atom 3
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at3, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), -0.6706364578997929, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), 0.5669584550606434, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), -0.7837441709115448, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), -0.004864516608601379, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), 0.02443134139590794, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), 0.02183605109402775, 1e-12 );
		}
		{
		//Atom 4
		Vector f1( 0.0 ), f2( 0.0 );
		dcst2.fill_f1_f2( at4, ConformationXYZ( ubqstump->conformation()), f1, f2, weights );
		TS_ASSERT_DELTA( f1.x(), 0.7154976591772094, 1e-12 );
		TS_ASSERT_DELTA( f1.y(), -0.6073138162523358, 1e-12 );
		TS_ASSERT_DELTA( f1.z(), 0.8278851575407408, 1e-12 );
		TS_ASSERT_DELTA( f2.x(), 0.005420850915159731, 1e-12 );
		TS_ASSERT_DELTA( f2.y(), -0.0263393747372001, 1e-12 );
		TS_ASSERT_DELTA( f2.z(), -0.02400679870736394, 1e-12 );
		}
		}

	}

	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void test_dihe_constraints()
	{
		using namespace core::scoring::constraints;
		using core::chemical::ResidueType;
		using core::chemical::AtomIndices;
		using core::conformation::Residue;
		using core::id::AtomID;

		pose::Pose start_pose( create_test_in_pdb_pose() );
		//core::import_pose::pose_from_pdb( start_pose, "core/scoring/constraints/test_in.pdb" );
		pose::Pose pose = start_pose; // a copy

		scoring::ScoreFunctionOP scorefxn = new scoring::ScoreFunction;
		scorefxn->reset();
		scorefxn->set_weight( scoring::fa_atr, 0.80 );
		scorefxn->set_weight( scoring::fa_rep, 0.44 );
		scorefxn->set_weight( scoring::fa_sol, 0.65 );
		scorefxn->set_weight( scoring::dihedral_constraint, 1.0 );

		// input stddev is in degrees, but dihedral constraints deal in radians
		core::Real const stddev_radians = numeric::conversions::radians( 5.0 );
		ConstraintSetOP constraints = new ConstraintSet();

		core::Size const which_res = 116; // Leu 164
		core::Size const which_chi = 2;
		Residue const & rsd = pose.residue(which_res);
		ResidueType const & rsd_type = pose.residue_type(which_res);
		TR << "Constraint on rsd " << rsd.name3() << " " << which_res << ", chi " << which_chi << std::endl;

		core::Real const start_chi_degrees = rsd.chi(which_chi);
		core::Real const start_chi_radians = numeric::conversions::radians( start_chi_degrees );
		FuncOP restr_func = new CircularHarmonicFunc( start_chi_radians, stddev_radians );
		AtomIndices chi_idx = rsd_type.chi_atoms(which_chi); // 1-based
		ConstraintOP constraint = new DihedralConstraint(
			AtomID(chi_idx[1], which_res),
			AtomID(chi_idx[2], which_res),
			AtomID(chi_idx[3], which_res),
			AtomID(chi_idx[4], which_res),
			restr_func
		);
		constraints->add_constraint( constraint );
		pose.constraint_set( constraints );
		TR << "Constraint: " << constraint->atom(1) << " " << constraint->atom(2) << " " << constraint->atom(3) << " " << constraint->atom(4) << std::endl;
		TR << "Starting position " << rsd.chi(which_chi) << " == " << start_chi_degrees << " degrees" << std::endl;

		(*scorefxn)( pose );
		TR << "Starting constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: 0)" << std::endl;
		TS_ASSERT_DELTA( 0.0, pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-6 );

		//rsd.chi( rsd.chi()+20.0 );
		pose.set_chi( which_chi, which_res, rsd.chi(which_chi)+20.0 );
		TR << "New position " << rsd.chi(which_chi) << " degrees" << std::endl;
		(*scorefxn)( pose );
		TR << "New constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: >0)" << std::endl;
		TS_ASSERT( 10.0 < pose.energies().total_energies()[ scoring::dihedral_constraint ] );

		kinematics::MoveMapOP mm = new kinematics::MoveMap;
		mm->set_chi( which_res, true );

		TR << "Minimizing..." << std::endl;
	// BAD!  We are in core, not protocols!  protocols lib is not even linked here!
	//protocols::moves::MinMover min_mover( mm, scorefxn, "dfpmin_armijo_nonmonotone_atol", 0.001, true /*use_nblist*/ );

		core::optimization::AtomTreeMinimizer minimizer;
		core::optimization::MinimizerOptionsOP min_options =
			new core::optimization::MinimizerOptions( "dfpmin_armijo_nonmonotone_atol", 0.001, true, true, false );
		minimizer.run( pose, *mm, *scorefxn, *min_options );

	//	min_mover.apply( pose );
		(*scorefxn)( pose );
		TR << "New position " << rsd.chi(which_chi) << " degrees" << std::endl;
		TS_ASSERT_DELTA( start_chi_degrees, rsd.chi(which_chi), 0.5 );
		TR << "New constraint energy " << pose.energies().total_energies()[ scoring::dihedral_constraint ] << " (expected: ~0)" << std::endl;
		TS_ASSERT_DELTA( 0.0, pose.energies().total_energies()[ scoring::dihedral_constraint ], 1e-1 );

		//[old behavior] When a pose is copied, we make a (very) deep copy of its constraints.
		// When a pose is copied, we make a shallow copy of its constraint set
		// Check that all constraints are copied and compare as equal.
		pose::Pose pose_copy( pose );
		TS_ASSERT( pose.constraint_set() != pose_copy.constraint_set() );
		utility::vector1< ConstraintCOP > constraints_orig = pose.constraint_set()->get_all_constraints();
		utility::vector1< ConstraintCOP > constraints_copy = pose_copy.constraint_set()->get_all_constraints();
		TS_ASSERT( constraints_orig.size() == 1 );
		TS_ASSERT( constraints_orig.size() == constraints_copy.size() );
		for(core::Size i = 1; i <= constraints_orig.size(); ++i) {
			TR << "Original " << i << ": " << constraints_orig[i]->to_string() << std::endl;
			TR << "Copy " << i << ": " << constraints_copy[i]->to_string() << std::endl;
			//[old behavior] Objects are equivalent (same data) but not identical (different places in memory)
			//TS_ASSERT( constraints_orig[i] != constraints_copy[i] );
			// Objects are identical (same place in memory)
			TS_ASSERT( constraints_orig[i] == constraints_copy[i] );
			TS_ASSERT( constraints_orig[i]() == constraints_copy[i]() );
			TS_ASSERT( constraints_orig[i]->to_string() == "DihedralConstraint 116 2 116 6 116 7 116 8 CircularHarmonicFunc 2.44839 0.0872665" );
			TS_ASSERT( constraints_orig[i]->to_string() == constraints_copy[i]->to_string() );
		}

		core::Size const size_before = pose.constraint_set()->get_all_constraints().size();
		TS_ASSERT( size_before == 1 );
		TS_ASSERT( pose.remove_constraint(constraint) == true );
		core::Size const size_after = pose.constraint_set()->get_all_constraints().size();
		TS_ASSERT( size_after == 0 );
		// Can't remove a constraint twice:
		TS_ASSERT( pose.remove_constraint(constraint) == false );
	}


};
