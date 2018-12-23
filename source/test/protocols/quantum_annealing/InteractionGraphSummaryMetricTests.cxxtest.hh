// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/quantum_annealing/InteractionGraphSummaryMetricTests.cxxtest.hh
/// @brief  Unit tests for the string metric that dumps out a summary of an interaction graph, for use by external packers.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/pdb1rpb.hh>

// Project Headers
#include <protocols/quantum_annealing/InteractionGraphSummaryMetric.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>

// Additional Protocols Headers
#include <protocols/quantum_annealing/ExternalPackerResultLoader.hh>


// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("InteractionGraphSummaryMetricTests");

#define CHI_DELTA 1e-2 //Error for chi-value comparisons

class InteractionGraphSummaryMetricTests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Set up a simple packing problem involving three residues from the trp cage (7*2*2=28 rotamer combinations), write it out, and
	/// read it back in.
	void test_summary_from_trpcage_repack() {
		using namespace protocols::quantum_annealing;

		core::pose::Pose pose( create_trpcage_ideal_pose() );

		core::pack::task::TaskFactoryOP task( utility::pointer::make_shared< core::pack::task::TaskFactory >() );
		//Prevent design:
		core::pack::task::operation::RestrictToRepackingOP nodesign( utility::pointer::make_shared< core::pack::task::operation::RestrictToRepacking >() );
		task->push_back( nodesign );
		//Prevent repacking except at residues 3, 6, and 7:
		core::select::residue_selector::ResidueIndexSelectorOP selection( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( "3,6,7" ) );
		core::pack::task::operation::PreventRepackingRLTOP restriction( utility::pointer::make_shared< core::pack::task::operation::PreventRepackingRLT >() );
		core::pack::task::operation::OperateOnResidueSubsetOP subset( utility::pointer::make_shared< core::pack::task::operation::OperateOnResidueSubset >( restriction, selection, true ) );
		task->push_back( subset );

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );

		InteractionGraphSummaryMetricOP summary( utility::pointer::make_shared< InteractionGraphSummaryMetric >() );
		summary->set_scorefunction( sfxn );
		summary->set_task_factory( task );
		std::string const problem_string( summary->calculate(pose) );

		TR << "\n" << problem_string << std::endl;

		//Reconstruct the solution with the mover:
		core::pose::Pose reconstructed_pose;
		std::string const solution_string( "3\t6\n6\t2\n7\t2" );
		protocols::quantum_annealing::ExternalPackerResultLoader result_loader;
		result_loader.set_packer_problem_string( problem_string );
		result_loader.set_rotamer_selection_string( solution_string );
		result_loader.apply( reconstructed_pose );

		TS_ASSERT_EQUALS( reconstructed_pose.total_residue(), 20 );
		//These expected chi values will have to be updated if the rotamer libraries change:
		TS_ASSERT_DELTA( reconstructed_pose.chi( 1, 3 ), -172.045, CHI_DELTA );
		TS_ASSERT_DELTA( reconstructed_pose.chi( 2, 3 ),  98.6754, CHI_DELTA );
		TS_ASSERT_DELTA( reconstructed_pose.chi( 1, 6 ), -177.923, CHI_DELTA );
		TS_ASSERT_DELTA( reconstructed_pose.chi( 2, 6 ),  62.6269, CHI_DELTA );
		TS_ASSERT_DELTA( reconstructed_pose.chi( 1, 7 ), -176.864, CHI_DELTA );
		TS_ASSERT_DELTA( reconstructed_pose.chi( 2, 7 ),  58.5823, CHI_DELTA );

		//reconstructed_pose.dump_pdb("QTEMP.pdb");
	}

};
