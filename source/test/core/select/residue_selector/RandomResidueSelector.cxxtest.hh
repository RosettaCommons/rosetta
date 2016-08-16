// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/RandomResidueSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::RandomResidueSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/RandomResidueSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers

using namespace core::select::residue_selector;


class RandomResidueSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_RandomResidueSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<RandomResidue name=allala num_residues=5 />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new RandomResidueSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueSubset subset = name_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.total_residue() );

		// test
		ResidueVector residues( subset );
		TS_ASSERT_EQUALS( residues.size(), 5 );

		// run a second time, and residues should be different
		ResidueVector residues2( name_rs->apply( trpcage ) );
		TS_ASSERT_EQUALS( residues2.size(), residues.size() );

		bool same = true;
		for ( ResidueVector::const_iterator r1=residues.begin(), r2=residues2.begin(); ( r1!=residues.end() ) && (r2!=residues2.end()); ++r1, ++r2 ) {
			if ( *r1 != *r2 ) {
				same = false;
				break;
			}
		}
		TS_ASSERT( !same );
	}

};
