// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/analysis/simple_metrics/SequenceRecoveryMetricTests.cxxtest.hh
/// @brief  Test the SequenceRecoveryMetrics
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/analysis/simple_metrics/SequenceRecoveryMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("SequenceRecoveryMetricTests");


class SequenceRecoveryMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_basic() {
		core::pose::PoseOP ref_pose( fullatom_poseop_from_string( trp_cage_ideal() ) );
		TR << "Sequence 1: " << ref_pose->sequence() << std::endl;

		core::pose::Pose pose( *ref_pose );
		TS_ASSERT( pose.residue(3).aa() != pose.residue(6).aa() );
		core::conformation::Residue tyr( pose.residue(3) );
		pose.replace_residue( 6, tyr, true );
		TR << "Sequence 2: " << pose.sequence() << std::endl;

		core::select::residue_selector::ResidueSelectorOP select( new core::select::residue_selector::TrueResidueSelector );

		protocols::analysis::simple_metrics::SequenceRecoveryMetric metric;

		metric.set_comparison_pose( ref_pose );
		metric.set_residue_selector( select );

		core::Real nres = pose.size();

		TS_ASSERT_EQUALS( metric.calculate( pose ),  (nres-1.0)/nres );

		TS_ASSERT( pose.residue(4).aa() != pose.residue(9).aa() );
		core::conformation::Residue ile( pose.residue(4) );
		pose.replace_residue( 9, ile, true );
		TR << "Sequence 3: " << pose.sequence() << std::endl;

		TS_ASSERT_EQUALS( metric.calculate( pose ),  (nres-2.0)/nres );
	}

	void test_subsection() {
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "NLYIQWLKDGGPSSGRPPPS", "fa_standard");

		core::select::residue_selector::ResidueSelectorOP select( new core::select::residue_selector::ResidueIndexSelector("4-15") ); // 12 residues inclusive

		protocols::analysis::simple_metrics::SequenceRecoveryMetric metric;
		metric.set_comparison_pose( core::pose::PoseOP( new core::pose::Pose(pose) ) );
		metric.set_residue_selector( select );

		TS_ASSERT_EQUALS( metric.calculate( pose ),  1.0 ); // No difference between them,

		//                                         .    *  * *   * . .
		core::pose::make_pose_from_sequence(pose, "ALYIQYLKNGAPSSTRQPSS", "fa_standard"); // 3 mutations outside range 4 inside.

		TS_ASSERT_EQUALS( metric.calculate( pose ), (12.0 - 4.0)/12.0 );
	}

	void test_pssm() {
		core::pose::Pose pose( fullatom_pose_from_string( trp_cage_ideal() ) );
		core::conformation::Residue leu( pose.residue(2) );
		core::conformation::Residue arg( pose.residue(16) );

		pose.replace_residue( 5, leu, true ); // These are negatives on the profile matrix
		pose.replace_residue(19, arg, true );

		pose.replace_residue( 4, leu, true ); // These are positives on the profile matrix - they shouldn't be counted as mutations
		pose.replace_residue( 8, arg, true );

		core::select::residue_selector::ResidueSelectorOP select( new core::select::residue_selector::TrueResidueSelector );

		protocols::analysis::simple_metrics::SequenceRecoveryMetric metric;
		metric.set_residue_selector( select );
		metric.load_pssm( "protocols/analysis/simple_metrics/trp_cage.pssm" );

		core::Real nres = pose.size();

		TS_ASSERT_EQUALS( metric.calculate( pose ),  (nres-2.0)/nres );
	}

	void test_average() {
		core::select::residue_selector::ResidueSelectorOP select( new core::select::residue_selector::TrueResidueSelector );

		protocols::analysis::simple_metrics::SequenceRecoveryMetric metric;
		metric.set_residue_selector( select );
		metric.load_pssm( "protocols/analysis/simple_metrics/trp_cage.pssm" );
		metric.set_use_ave_pssm( true );

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "NLYIQWLKDGGPSSGRPPPS", "fa_standard");

		core::Real native( 3+3+7+3+5+10+4+4+5+6+6+7+4+4+6+5+7+7+7+4 );
		// Native
		TS_ASSERT_EQUALS( metric.calculate( pose ), native/pose.size() );

		core::pose::make_pose_from_sequence(pose, "NLYRQWLKDGGPSSGRPTPS", "fa_standard");

		TS_ASSERT_EQUALS( metric.calculate( pose ), (native - 3 + -2 - 7 + -1)/pose.size()  );
	}

	void test_delta() {
		core::select::residue_selector::ResidueSelectorOP select( new core::select::residue_selector::TrueResidueSelector );

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose, "NLYIQWLKDGGPSSGRPPPS", "fa_standard");

		protocols::analysis::simple_metrics::SequenceRecoveryMetric metric;
		metric.set_residue_selector( select );
		metric.load_pssm( "protocols/analysis/simple_metrics/trp_cage.pssm" );
		metric.set_use_ave_pssm( true );
		metric.set_comparison_pose( core::pose::PoseOP( new core::pose::Pose(pose) ) );

		// Native
		TS_ASSERT_EQUALS( metric.calculate( pose ), 0 );

		core::pose::make_pose_from_sequence(pose, "NLYRQWLKDGGPSSGRPTPW", "fa_standard");

		TS_ASSERT_EQUALS( metric.calculate( pose ), ((-2.0 - 3.0) + (-1.0 - 7.0) + (-3.0 - 4.0) )/pose.size()  );
	}



};
