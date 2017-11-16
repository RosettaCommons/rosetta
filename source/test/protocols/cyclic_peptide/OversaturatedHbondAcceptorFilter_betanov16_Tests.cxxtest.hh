// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter_betanov16_Tests.cxxtest.hh
/// @brief  Unit tests for a filter that flags excessive numbers of hydrogen bonds to acceptors.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.cyclic_peptide.OversaturatedHbondAcceptorFilter_betanov16_Tests");


class OversaturatedHbondAcceptorFilter_betanov16_Tests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init_with_additional_options("-no_optH -beta_nov16");

		oversaturated_bb_bb_bb_ = core::import_pose::pose_from_file( "protocols/cyclic_peptide/oversaturated_bb_bb_bb.pdb", false, core::import_pose::PDB_file);

		//Link the termini:
		protocols::cyclic_peptide::DeclareBond decbond1;
		decbond1.set( 1, "N", 8, "C", false );
		decbond1.apply(*oversaturated_bb_bb_bb_);

	}

	void tearDown(){
	}

	/// @brief Check that the filter flags an extra backbone-backbone hydrogen bond.
	///
	void test_bb_bb_extra_hbond() {
		using namespace protocols::cyclic_peptide;
		OversaturatedHbondAcceptorFilter filter;  filter.set_max_allowed_oversaturated(0);
		filter.set_consider_mainchain_only(false);

		TR << "Testing backbone-backbone oversaturation in all-hbonds mode." << std::endl;
		TS_ASSERT( !filter.apply( *oversaturated_bb_bb_bb_ ) );
		TS_ASSERT_DELTA( filter.report_sm(*oversaturated_bb_bb_bb_), 1, 0.00001 );

		filter.set_consider_mainchain_only(true);

		TR << "Testing backbone-backbone oversaturation in backbone-only mode." << std::endl;
		TS_ASSERT( !filter.apply( *oversaturated_bb_bb_bb_ ) );
		TS_ASSERT_DELTA( filter.report_sm(*oversaturated_bb_bb_bb_), 1, 0.00001 );
	}

private:

	core::pose::PoseOP oversaturated_sc_bb_bb_;

	core::pose::PoseOP oversaturated_bb_bb_bb_;

};



