// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreTest.cxxtest.hh
/// @brief  unified scoring test.
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/EnvPairPotential.hh>

// Project Headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

//#include <core/chemical/ResidueTypeSet.hh>


#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/types.hh>

// Package headers

#include <test/UTracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <utility/down_cast.hh>
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
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <execinfo.h>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.EnvPairPotential.cxxtest");

// using declarations
using namespace core;
using namespace scoring;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreTest
/// @brief: unified tests for difference score functions/methods
///////////////////////////////////////////////////////////////////////////
class EnvPairPotentialTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-no_optH" );
	}

	void tearDown() {}

	void test_compute_centroid_environment() {
		//using core::pose::datacache::CacheableDataType::CEN_LIST_INFO;
		pose::Pose pose;
		core::import_pose::centroid_pose_from_pdb( pose,"core/scoring/test_in.pdb", core::import_pose::PDB_file);
		Real fcen_gold[][3] = {
			{ 4, 16, 18 },
			{ 4, 10, 11 },
			{ 2, 8, 8 },
			{ 2.06444, 15, 17.9356 },
			{ 5, 18.6177, 25 },
			{ 2, 13, 13 },
			{ 3, 14, 17 },
			{ 5, 20, 33 },
			{ 2.91155, 17, 29.3436 },
			{ 2, 10, 14 },
			{ 2, 12, 20.2338 },
			{ 4, 22.1345, 28 },
			{ 3, 13, 21.3865 },
			{ 2, 8, 15 },
			{ 3, 16, 21.0463 },
			{ 6.22525, 18, 29.3577 },
			{ 2, 12, 19 },
			{ 2, 11, 16 },
			{ 4, 15, 18 },
			{ 3, 10.259, 20 },
			{ 6.42046, 19, 25.7885 },
			{ 2, 13, 20 },
			{ 4, 15, 19 },
			{ 3, 8, 9 },
			{ 2, 14, 21 },
			{ 3.34125, 21, 27.6587 },
			{ 3, 17, 17 },
			{ 2, 8, 13 },
			{ 4, 13, 20 },
			{ 4, 14, 15 },
			{ 4, 10, 14.8917 },
			{ 3, 10, 12 },
			{ 4, 13, 18 },
			{ 1, 8, 11 },
			{ 4, 11, 15 },
			{ 3.43142, 19.6354, 28.2873 },
			{ 3, 15, 24.9727 },
			{ 2, 10, 13 },
			{ 4.06444, 15, 19.9356 },
			{ 4.58881, 23, 32.4112 },
			{ 2, 13, 18 },
			{ 1, 8, 12 },
			{ 5, 17, 19.6263 },
			{ 5.06895, 17.9933, 22.0363 },
			{ 3, 10, 15.2058 },
			{ 4, 15, 17.92 },
			{ 4, 15, 13 },
			{ 5, 22, 34.8917 },
			{ 5, 18, 23.2058 },
			{ 2, 11.6297, 17 },
			{ 3.78291, 17.6177, 26.4783 },
			{ 5.92079, 21, 32.0792 },
			{ 2, 15, 26 },
			{ 2.78291, 17, 29.2171 },
			{ 4, 18.9997, 30.6509 },
			{ 3.20133, 25, 31.7987 },
			{ 4, 19.6297, 30.1747 },
			{ 3.73607, 19, 27.2639 },
			{ 4, 17, 31 },
			{ 5, 26, 31.6457 },
			{ 3.98838, 25, 30.2966 },
			{ 3, 12.0705, 19 },
			{ 4, 15.8523, 22 },
			{ 4, 21.5688, 31.926 },
			{ 1, 14, 24 },
			{ 2, 6, 13 },
			{ 3, 11, 15 },
			{ 1, 6, 13 },
			{ 7, 13, 14 },
			{ 3, 14.8933, 14.4096 },
			{ 7, 19, 28 },
			{ 2.66277, 12.8869, 16.3372 },
			{ 4, 10.8523, 12.6457 },
			{ 6, 21, 28 },
			{ 4.66277, 22.7379, 25.3372 },
			{ 1, 11.6884, 18 },
			{ 3, 14, 18.6562 },
			{ 4.22525, 22, 31.7748 },
			{ 3.20133, 17.1515, 22.7987 },
			{ 6, 21, 25 },
			{ 2, 7.79509, 10 },
			{ 3, 12, 14 },
			{ 4, 16, 16 },
			{ 2, 7.7379, 13 },
			{ 4, 12.259, 20.1747 },
			{ 5, 20, 34.1504 },
			{ 2.68237, 17, 24.3176 },
			{ 1, 8, 12 },
			{ 4, 17.1515, 26 },
			{ 5, 24, 31 },
			{ 2.68237, 13, 17.5875 },
			{ 1, 9, 18.2233 },
			{ 2, 16.2472, 28 },
			{ 5, 20.0949, 26.7475 },
			{ 1, 10.2472, 13 },
			{ 6, 15, 20 },
			{ 5, 14.8632, 16.2849 },
			{ 5, 22, 28 },
			{ 3, 12.8869, 18 },
			{ 3, 10, 14 },
			{ 4.25231, 23, 26.9644 },
			{ 4, 17, 26 },
			{ 1, 6.57749, 14 },
			{ 2, 11.0949, 16 },
			{ 1, 6.87186, 6 },
			{ 1, 5, 7 },
			{ 4, 15, 17.706 },
			{ 2, 12.0283, 15 },
			{ 2, 9, 12.706 },
			{ 2, 11, 13.3322 },
			{ 2, 13, 17 },
			{ 5, 15, 18 },
			{ 3, 13, 24.3322 },
			{ 2, 11, 12.8348 },
			{ 3, 10.0422, 16 },
			{ 2, 12.5055, 18 }};
		ScoreFunction sfxn;
		sfxn.set_weight( vdw, 1.0 );
		sfxn(pose);

		EnvPairPotential envp;
		envp.compute_centroid_environment( pose );
		TS_ASSERT( pose.data().has( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ) );
		CenListInfo const & cenlist( *( utility::pointer::static_pointer_cast< core::scoring::CenListInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ) )));

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			TS_ASSERT_DELTA( cenlist.fcen6( ii ),  fcen_gold[ ii -1 ][ 0 ], 1E-4 );
			TS_ASSERT_DELTA( cenlist.fcen10( ii ), fcen_gold[ ii -1 ][ 1 ], 1E-4 );
			TS_ASSERT_DELTA( cenlist.fcen12( ii ), fcen_gold[ ii -1 ][ 2 ], 1E-4 );
		}
	}

	void
	test_evaluate_cbeta_scores()
	{
		//using core::pose::datacache::CacheableDataType::CEN_LIST_INFO;

		Real cbeta_scores_gold[] = {
			0.280195, 0.294517, 0.581433, 0.53862, 0.244031,
			0.55759, 0.342203, 0.260113, 0.373973, 0.55647,
			0.55171, 0.289956, 0.342943, 0.555109, 0.342644,
			0.295597, 0.550949, 0.554389,0.280195, 0.341403,
			0.293721, 0.551429, 0.278995, 0.366926, 0.552629,
			0.330478, 0.342203, 0.55759, 0.279475, 0.283155,
			0.283303, 0.352684, 0.280195, 0.565751, 0.283155,
			0.325995, 0.346804, 0.55759, 0.276809, 0.275935,
			0.552149, 0.56199, 0.23841, 0.243273, 0.344935,
			0.280201, 0.285636, 0.265903, 0.241637, 0.552229,
			0.301598, 0.292971, 0.55919, 0.400759, 0.296656,
			0.348455, 0.295208, 0.305712, 0.297717, 0.257865,
			0.296298, 0.340923, 0.281555, 0.299199, 0.554709,
			0.55759, 0.345083, 0.55687, 0.305958, 0.345886,
			0.311399, 0.414463, 0.28745, 0.286436, 0.258601,
			0.551429, 0.341335, 0.289748, 0.33159, 0.281395,
			0.571431, 0.346443, 0.282435, 0.55759, 0.279685,
			0.265014, 0.412572, 0.56199, 0.287236, 0.256832,
			0.408867, 0.551161, 0.56191, 0.247786, 0.55687,
			0.275955, 0.240935, 0.249071, 0.342123, 0.346443,
			0.278772, 0.287236, 0.555749, 0.554389, 0.599915,
			0.587753, 0.280219, 0.555109, 0.559095, 0.557218,
			0.552229, 0.23931, 0.345882, 0.558436, 0.344363,
			0.552149 };
		ScoreFunction sfxn;
		sfxn.set_weight( vdw, 1.0 );
		pose::Pose pose;
		core::import_pose::centroid_pose_from_pdb( pose,"core/scoring/test_in.pdb");

		sfxn( pose );

		EnvPairPotential envp;
		envp.compute_centroid_environment( pose );
		TS_ASSERT( pose.data().has( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ) );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			Real env_e( 0.0 ), cb_score6( 0.0 ), cb_score12( 0.0 );
			envp.evaluate_env_and_cbeta_scores( pose, pose.residue( ii ),
				env_e, cb_score6, cb_score12 );
			Real cbscore = 2.667 * ( cb_score6 + cb_score12 ) * 0.3;
			TS_ASSERT_DELTA( cbscore, cbeta_scores_gold[ ii - 1 ], 1E-6 );

		}
	}


};
