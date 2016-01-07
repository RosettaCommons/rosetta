// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/loop_graph//LoopGraphTest.cxxtest.hh
/// @brief  testing LoopGraph class
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/scoring/loop_graph/LoopGraph.hh>
#include <core/scoring/loop_graph/LoopCycle.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>

#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("LoopGraphTest");


class LoopGraphTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( " -s core/scoring/loop_graph/GC_pair.pdb core/scoring/loop_graph/P1.pdb core/scoring/loop_graph/P2_part1.pdb core/scoring/loop_graph/P2_part2.pdb core/scoring/loop_graph/P3.pdb core/scoring/loop_graph/SRL_1JBT.pdb -fasta core/scoring/loop_graph/bound_seq_GC_SRL.fasta" );
	}

	void tearDown(){
	}

	void test_loop_cycles(){

		using namespace core;
		using namespace core::chemical;
		using namespace core::scoring::loop_graph;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
    using namespace protocols::stepwise::setup;
    using namespace utility::tools;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose = initialize_pose_and_other_poses_from_command_line( rsd_set );

		LoopGraph loop_graph;
		loop_graph.set_error_out_on_complex_cycles( false );
		loop_graph.update( *pose );
		utility::vector1< LoopCycle > const & loop_cycles = loop_graph.loop_cycles();
		TS_ASSERT_EQUALS( loop_cycles.size(), 6 );
		for ( Size n = 1; n <= loop_cycles.size(); n++ ) {
			TR << "#"<< n << ": " << loop_cycles[ n ] << std::endl;
		}

		// Loop definition is: takeoff_pos, landing_pos, takeoff_domain, landing_domain;
		TS_ASSERT( loop_cycles.has_value( LoopCycle( make_vector1( Loop( 11,16,4,4) ) ) ) );
		// following two loop cycles are *nested*, which is what is called a complex_loop_graph.
		TS_ASSERT( loop_cycles.has_value( LoopCycle( make_vector1( Loop(   23, 25, 1, 6 ),
																															 Loop(   53, 55, 6, 2 ),
																															 Loop(   5, 7, 2, 3 ),
																															 Loop(  20, 23, 3, 1 ) ) ) ) );
		TS_ASSERT( loop_cycles.has_value( LoopCycle( make_vector1( Loop(   23, 25, 1, 6 ),
																															 Loop(   53, 55, 6, 2 ),
																															 Loop(   59, 60, 2, 1 ) ) ) ) );
		TS_ASSERT( !loop_graph.has_just_simple_cycles() );




	}



};



