// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/residue_selectors/PairedSheetResidueSelectorTests.cxxtest.hh
/// @brief  Test suite for PairedSheetResidueSelector
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.hh>


// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR("PairedSheetResidueSelectorTests");

class PairedSheetResidueSelectorTests : public CxxTest::TestSuite {
	//Define Variables
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;
	typedef std::vector< core::Size > ResidueVector;
	typedef protocols::denovo_design::residue_selectors::PairedSheetResidueSelector PairedSheetResidueSelector;

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief checks whether the residues in the given vector are true in the subset.  All other residues should be false.
	void
	check_true( ResidueSubset const & subset, ResidueVector const & true_residues ) const
	{
		ResidueVector::const_iterator next_paired_res = true_residues.begin();
		ResidueVector::const_iterator end = true_residues.end();
		for ( core::Size resid=1; resid<=subset.size(); ++resid ) {
			if ( (next_paired_res!=end) && ( *next_paired_res == resid ) ) {
				TR.Debug << "Found expected paired res: " << resid << std::endl;
				TS_ASSERT( subset[ resid ] );
				++next_paired_res;
			} else {
				TS_ASSERT( !subset[ resid ] );
			}
		}
	}

	void test_paired_res()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/fldsgn/filters/test_sheet.pdb" );

		PairedSheetResidueSelector selector;
		selector.set_sheet_topology( "1-2.P.1" );

		ResidueSubset subset = selector.apply( pose );
		TS_ASSERT_EQUALS( subset.size(), pose.total_residue() );

		// strand 1 is 2-6, but only 3-6 should be paired
		// strand 2 is 32-38, but only 32-35 should be paired
		ResidueVector should_be_paired = boost::assign::list_of (3)(4)(5)(6)(32)(33)(34)(35);
		check_true( subset, should_be_paired );

		// changing the register shift should change the paired residues
		selector.set_sheet_topology( "1-2.P.-1" );
		subset = selector.apply( pose );
		TS_ASSERT_EQUALS( subset.size(), pose.total_residue() );

		// strand 1 is 2-6, and 2-6 should be paired
		// strand 2 is 32-38, and 33-37 should be paired
		should_be_paired = boost::assign::list_of (2)(3)(4)(5)(6)(33)(34)(35)(36)(37).convert_to_container< ResidueVector >();
		check_true( subset, should_be_paired );

		// new secondary structure adds a strand residue to E1
		std::string secstruct = "LEEEEEEHHHHLLLLHHHHHHHHHHHHLLLLEEEEEEELLLHHHHHHHHHHHHHHHHHHLLLEEEEEEEL";
		selector.set_secstruct( secstruct );
		subset = selector.apply( pose );
		TS_ASSERT_EQUALS( subset.size(), pose.total_residue() );

		// strand 1 is 2-7, and 2-7 should be paired
		// strand 2 is 32-38, and 33-38 should be paired
		should_be_paired = boost::assign::list_of (2)(3)(4)(5)(6)(7)(33)(34)(35)(36)(37)(38).convert_to_container< ResidueVector >();
		check_true( subset, should_be_paired );

		/// try adding both strands
		selector.set_sheet_topology( "1-2.P.1;2-3.P.0" );
		selector.set_secstruct( "" );
		subset = selector.apply( pose );
		TS_ASSERT_EQUALS( subset.size(), pose.total_residue() );

		// strand 1 is 2-6, and 3-6 should be paired
		// strand 2 is 32-38, and 32-38 should be paired
		// strand 3 is 63-69 and all should be paired
		should_be_paired = boost::assign::list_of (3)(4)(5)(6)(32)(33)(34)(35)(36)(37)(38)(63)(64)(65)(66)(67)(68)(69).convert_to_container< ResidueVector >();
		check_true( subset, should_be_paired );
	}

};



