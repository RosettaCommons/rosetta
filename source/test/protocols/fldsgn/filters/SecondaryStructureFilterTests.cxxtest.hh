// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/fldsgn/filters/SecondaryStructureFilterTests.cxxtest.hh
/// @brief  Test suite for Secondary structure filter
/// @author Tom Linsky (tlinsky@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SecondaryStructureFilterTests");

namespace ss_test {
class SecondaryStructureFilterProt : public protocols::fldsgn::filters::SecondaryStructureFilter {

public:
	void
	modify_secstruct( Pose const & pose, String & secstruct ) const
	{ correct_for_incomplete_strand_pairings( pose, secstruct ); }

	std::string
	get_filtered_secstruct_( Pose const & pose ) const
	{ return get_filtered_secstruct( pose ); }

	core::Size
	compute_( Pose const & pose, core::select::residue_selector::ResidueSubset const & subset ) const
	{ return compute( pose, subset ); }

};
}

class SecondaryStructureFilterTests : public CxxTest::TestSuite {
public:
	typedef ss_test::SecondaryStructureFilterProt SecondaryStructureFilter;
	//Define Variables

public:

	void setUp()
	{
		core_init();
	}

	void tearDown()
	{
	}

	void test_convert_impossible_pairings()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/fldsgn/filters/ss_filter_test.pdb" );

		// set pose ss
		std::string const built_ss =            "LEEEEEELLLHHHHHHHHHHHHHHHLLEEEEEL";
		std::string const bad_ss =              "LEEEEEELLLHHHHHHHHHHHHHLLLLEEEEEL";
		std::string const modified_ss =         "LhEEEEELLLHHHHHHHHHHHHHHHLLEEEEEL";
		std::string const modified_ss_shifted = "LEEEEhhLLLHHHHHHHHHHHHHHHLLhEEEEL";

		core::Size resid = 1;
		for ( std::string::const_iterator ss=built_ss.begin(); ss!=built_ss.end(); ++ss, ++resid ) {
			pose.set_secstruct( resid, *ss );
		}

		core::scoring::dssp::Dssp dssp( pose );
		TR << dssp.get_dssp_secstruct() << " " << pose.size() << std::endl;

		SecondaryStructureFilter filt;
		filt.set_use_dssp( true );
		filt.filtered_ss( built_ss );

		// filtered ss should be correctly determined
		TS_ASSERT_EQUALS( filt.get_filtered_secstruct_( pose ), built_ss );

		// we should fail without setting strand pairings - DSSP will never recognize residue 2 as E
		TS_ASSERT( !filt.apply( pose ) );

		/// try to add a strand pairing
		std::string const spairstr = "1-2.P.1";
		filt.set_strand_pairings( spairstr );

		// built ss should be equal to get_filtered_secstruct_
		std::string filt_ss = filt.get_filtered_secstruct_( pose );
		TS_ASSERT_EQUALS( filt_ss, built_ss );

		// modified ss should be equal to built_ss since there are no strand pairings
		filt.modify_secstruct( pose, filt_ss );
		TS_ASSERT_EQUALS( filt_ss, modified_ss );

		// filter should now pass
		TS_ASSERT( filt.apply( pose ) );

		// changing wanted secondary structure should cause failure
		filt.set_use_dssp( false );
		std::string::const_iterator ss = bad_ss.begin();
		for ( core::Size resid=1; resid<=pose.size(); ++resid, ++ss ) {
			pose.set_secstruct( resid, *ss );
		}
		TS_ASSERT_EQUALS( filt.report_sm(pose),  core::Real(31)/core::Real(33) );
		filt.set_use_dssp( true );

		/// try to add a different strand pairing
		filt.set_strand_pairings( "1-2.P.-1" );

		filt_ss = filt.get_filtered_secstruct_( pose );
		TS_ASSERT_EQUALS( filt_ss, built_ss );

		filt.modify_secstruct( pose, filt_ss );
		TS_ASSERT_EQUALS( filt_ss, modified_ss_shifted );

		// filter should fail because residue 2 is supposed to be paired and E but isn't
		TS_ASSERT( !filt.apply( pose ) );
	}

};
