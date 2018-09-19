// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/core/simple_metrics/metrics/SelectedResidueCountMetricTests.cxxtest.hh
/// @brief  Unit tests for the SelectedResidueCountMetric.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/metrics/SelectedResidueCountMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Protocols Headers

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("SelectedResidueCountMetricTests");


class SelectedResidueCountMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}


	/// @brief Confirm that the counter counts the number of residues in the selection.
	void test_with_selector(){
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/simple_metrics/metrics/testpose.pdb", false, core::import_pose::PDB_file );

		core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		index_selector->set_index("3,5,7-12"); //Eight residues selected

		core::simple_metrics::metrics::SelectedResidueCountMetric selcount;
		selcount.set_residue_selector(index_selector);
		core::Real const count( selcount.calculate(pose) );
		TS_ASSERT_EQUALS(count, 8);
	}

	/// @brief Confirm that the counter counts the number of residues in the pose, in the absence of a selector.
	void test_without_selector(){
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/simple_metrics/metrics/testpose.pdb", false, core::import_pose::PDB_file );
		core::simple_metrics::metrics::SelectedResidueCountMetric selcount;
		core::Real const count( selcount.calculate(pose) );
		TS_ASSERT_EQUALS(count, pose.total_residue());
	}



};
