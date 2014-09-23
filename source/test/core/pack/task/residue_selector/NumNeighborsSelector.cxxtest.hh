// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/NumNeighborsSelector.cxxtest.hh
/// @brief  test suite for core::pack::task::residue_selector::NumNeighborsSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/pack/task/residue_selector/DummySelectors.hh>

// Package headers
#include <core/pack/task/residue_selector/NumNeighborsSelector.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task::residue_selector;


class NumNeighborsSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief test that we select all the residues in trpcage
	void test_num_neighbors_selector_1() {
		NumNeighborsSelectorOP nneighbs_rs = new NumNeighborsSelector( 12, 10.0 );
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset( trpcage.total_residue(), false );
		nneighbs_rs->apply( trpcage, subset );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			core::conformation::Residue const & iires = trpcage.residue(ii);
			core::Size count_neighbs = 0;
			for ( core::Size jj = 1; jj <= trpcage.total_residue(); ++jj ) {
				if ( ii == jj ) continue;
				core::conformation::Residue const & jjres = trpcage.residue(jj);
				if ( iires.xyz( iires.nbr_atom() ).distance_squared( jjres.xyz( jjres.nbr_atom() ) ) <= 100.0 ) ++count_neighbs;
			}
			TS_ASSERT_EQUALS( subset[ ii ], count_neighbs >= 12 );
		}
	}


	/// @brief Test NumNeighborsSelector::parse_my_tag using default options
	void test_num_neighbors_selector_parse_my_tag_w_defaults() {
		std::string tag_string = "<NumNeighbors name=nneighbs_rs/>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		NumNeighborsSelectorOP nn_rs = new NumNeighborsSelector;
		try {
			nn_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
			return;
		}
		TS_ASSERT_EQUALS( nn_rs->count_water(), false );
		TS_ASSERT_EQUALS( nn_rs->threshold(), core::Size(17) );
		TS_ASSERT_EQUALS( nn_rs->distance_cutoff(), 10.0 );
	}

	/// @brief Test NumNeighborsSelector::parse_my_tag with provided options
	void test_num_neighbors_selector_parse_my_tag() {
		std::string tag_string = "<NumNeighbors name=nneighbs_rs count_water=false threshold=8 distance_cutoff=5.5 />";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		NumNeighborsSelectorOP nn_rs = new NumNeighborsSelector;
		try {
			nn_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
			return;
		}
		TS_ASSERT_EQUALS( nn_rs->count_water(), false );
		TS_ASSERT_EQUALS( nn_rs->threshold(), core::Size(8) );
		TS_ASSERT_EQUALS( nn_rs->distance_cutoff(), 5.5 );
	}

};
