// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methodsDEEREnergyMethod_test.cxxtest.hh
/// @brief  Tests for DEEREnergyMethod class
/// @author Diego del Alamo

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/energy_methods/DEEREnergy.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <numeric/xyzVector.hh>

#include <map>

class DEEREnergyMethodTest : public CxxTest::TestSuite {
public:

	void setUp() {
		std::string data_file = "-epr_deer:input_files core/scoring/epr_deer/input_file_epr_deer.txt";
		core_init_with_additional_options( data_file );
	}

	void tearDown() {
		// NULL
	}

	void test_energy_method_big() {
		core::pose::Pose test_pose;
		core::io::pdb::build_pose_from_pdb_as_is( test_pose, "core/scoring/epr_deer/2lzm.pdb" );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight_if_zero( core::scoring::ScoreType::epr_deer_score, 1.0 );
		sfxn.setup_for_scoring( test_pose );
		core::scoring::epr_deer::DEERDataCacheOP datacache = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDataCache >(
			test_pose.data().get_ptr( core::pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
		TS_ASSERT_EQUALS( datacache->size(), 6 );
		TS_ASSERT_EQUALS( datacache->edge( 60, 94 ).size(), 5 );
		TS_ASSERT_EQUALS( datacache->edge( 94, 123 ).size(), 1 );
	}

	void test_energy_method_details() {
		core::pose::Pose test_pose;
		core::io::pdb::build_pose_from_pdb_as_is( test_pose, "core/scoring/epr_deer/2lzm.pdb" );

		core::energy_methods::DEEREnergy energy_method;
		energy_method.initialize_energy_method( test_pose );
		core::scoring::epr_deer::DEERDataCacheOP datacache = utility::pointer::dynamic_pointer_cast< core::scoring::epr_deer::DEERDataCache >(
			test_pose.data().get_ptr( core::pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
		core::scoring::epr_deer::EPRSpinLabel sl;
		core::Real score_1 = energy_method.get_score( test_pose, sl, 1 );
		TS_ASSERT_EQUALS( score_1, ( *datacache )[ 1 ]->score() );
		TS_ASSERT_EQUALS( energy_method.defines_residue_pair_energy( test_pose, 94, 123 ), true );
		TS_ASSERT_EQUALS( energy_method.defines_residue_pair_energy( test_pose, 94, 155 ), false );

		core::scoring::ScoreFunction sfxn;
		energy_method.setup_for_derivatives( test_pose, sfxn );
		TS_ASSERT( datacache->f1_force( 60 ) != numeric::xyzVector< core::Real >( 0.0, 0.0, 0.0 ) );
		TS_ASSERT( datacache->f2_force( 60 ) != numeric::xyzVector< core::Real >( 0.0, 0.0, 0.0 ) );
	}
};
