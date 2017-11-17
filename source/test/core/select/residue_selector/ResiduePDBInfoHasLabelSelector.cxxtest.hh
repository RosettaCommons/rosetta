// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/ResiduePDBInfoHasLabelSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::ResiduePDBInfoHasLabelSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers
#include <string>

using namespace core::select::residue_selector;


class ResiduePDBInfoHasLabelSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test NotResidueSelector::parse_my_tag
	void test_ResiduePDBInfoHasLabelSelector_parse_my_tag() {
		std::stringstream ss;
		ss << "<ResiduePDBInfoHasLabel name=allala property=\"TESTPROPERTY\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResiduePDBInfoHasLabelSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
		} catch (utility::excn::Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		trpcage.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( trpcage ) ) );
		std::set< core::Size > residues = boost::assign::list_of (9)(20)(2)(7);
		for ( std::set< core::Size >::const_iterator res = residues.begin(); res != residues.end(); ++res ) {
			trpcage.pdb_info()->add_reslabel( *res, "TESTPROPERTY" );
		}

		ResidueSubset subset = name_rs->apply( trpcage );
		TS_ASSERT_EQUALS( subset.size(), trpcage.size() );

		// test
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			TS_ASSERT( !subset[ ii ] || ( residues.find( ii ) != residues.end() ) );
		}
	}

	// make sure we fail if no selection string is provided
	void test_NeighbohoodResidueSelector_fail_no_label() {
		std::stringstream ss;
		ss << "<ResiduePDBInfoHasLabel name=I_should_fail />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		ResidueSelectorOP name_rs( new ResiduePDBInfoHasLabelSelector );
		try {
			name_rs->parse_my_tag( tag, dm );
			TS_ASSERT( false ); //parsing should fail!
		} catch (utility::excn::Exception e ) {
			// should fail with a message
			TS_ASSERT( e.msg() != "" );
		}
	}

};
