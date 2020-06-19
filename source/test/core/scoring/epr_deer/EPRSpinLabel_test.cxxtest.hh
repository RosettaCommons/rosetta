// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/EPRSpinLabel_test.cxxtest.hh
/// @brief  Tests for EPRSpinLabel class
/// @author Diego del Alamo

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <map>

class EPRSpinLabelTest : public CxxTest::TestSuite {

public:

	void setUp() {
		std::string data_file = "-epr_deer:input_files core/scoring/epr_deer/input_file_epr_deer.txt";
		std::string coords_file = "-epr_deer:coords_files core/scoring/epr_deer/coords_file.txt";
		core_init_with_additional_options( data_file + " " + coords_file );
	}

	void tearDown() {
		// NULL
	}

	void test_histograms() {
		core::pose::Pose test_pose;
		core::io::pdb::build_pose_from_pdb_as_is( test_pose, "core/scoring/epr_deer/2lzm.pdb" );
		test_pose.update_residue_neighbors();

		core::scoring::epr_deer::EPRSpinLabel sl;
		std::pair< core::Size, std::string > res1 = std::make_pair( 60, "DEFAULT" );
		std::pair< core::Size, std::string > res2 = std::make_pair( 94, "DEFAULT_FAST" );
		for ( auto & item : { res1, res2 } ) {
			sl.label( item.first, item.second, test_pose );
		}
		auto histogram_1 = sl.histogram( { res1, res2 }, 2 );
		TS_ASSERT( histogram_1.size() > 0 );
		auto histogram_2 = sl.histogram( res1, res2, 2 );
		TS_ASSERT( histogram_2.size() > 0 );

		// spot check
		TS_ASSERT_DELTA( sl.normalize_distribution( histogram_1 )[ 49 ], sl.normalize_distribution( histogram_2 )[ 49 ], 1e-6 );

		sl.label( 123, "DEFAULT", test_pose );
		TS_ASSERT_EQUALS( sl[ std::make_pair( 123, "DEFAULT" ) ].size(), 50 );
	}

};
