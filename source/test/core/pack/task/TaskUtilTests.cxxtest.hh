// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/task/UtilTests.cxxtest.hh
/// @brief  Test utility functions for task utils.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("UtilTests");

using namespace core::pack::task;

class TaskUtilTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_res_agnostic_command_parsing() {

		///Check to make sure all the parsing actually works.
		//  Tests for TaskOp will make sure the classes are actually initialized properly.
		std::map< std::string, ResfileCommandOP > cmd_map = create_command_map();

		std::string cmd_str = "AROMATIC";
		utility::vector1< ResfileCommandOP > commands = parse_res_agnostic_commands( cmd_str, cmd_map );

		TS_ASSERT( commands.size()  == 1 );
		TS_ASSERT_EQUALS( commands[1]->get_name(), "AROMATIC" );


		cmd_str = "PIKAA ST";
		commands = parse_res_agnostic_commands( cmd_str, cmd_map );

		TS_ASSERT( commands.size()  == 1 );
		TS_ASSERT_EQUALS( commands[1]->get_name(), "PIKAA" );

	}

	void test_sequence_motif_parsing() {

		std::string motif = "N[^P][ST]-X";

		utility::vector1< utility::vector1< ResfileCommandOP > > commands = get_resfile_commands( motif );

		TS_ASSERT_EQUALS( commands.size(), 5 );

		//Make sure each position only has one command.
		TS_ASSERT_EQUALS( commands[1].size(), 1 );
		TS_ASSERT_EQUALS( commands[2].size(), 1 );
		TS_ASSERT_EQUALS( commands[3].size(), 1 );
		TS_ASSERT_EQUALS( commands[4].size(), 1 );
		TS_ASSERT_EQUALS( commands[5].size(), 1 );

		//Make sure each command is the correct type.
		TS_ASSERT_EQUALS( commands[1][1]->get_name(), "PIKAA" );
		TS_ASSERT_EQUALS( commands[2][1]->get_name(), "NOTAA" );
		TS_ASSERT_EQUALS( commands[3][1]->get_name(), "PIKAA" );
		TS_ASSERT_EQUALS( commands[4][1]->get_name(), "NATAA" );
		TS_ASSERT_EQUALS( commands[5][1]->get_name(), "ALLAA" );

	}

	void test_linear_sequence_motif_parsing() {

		std::string motif = "N-T";

		utility::vector1< utility::vector1< ResfileCommandOP > > commands = get_resfile_commands( motif );

		TS_ASSERT_EQUALS( commands.size(), 3 );

		//Make sure each position only has one command.
		TS_ASSERT_EQUALS( commands[1].size(), 1 );
		TS_ASSERT_EQUALS( commands[2].size(), 1 );
		TS_ASSERT_EQUALS( commands[3].size(), 1 );

		//Make sure each command is the correct type.
		TS_ASSERT_EQUALS( commands[1][1]->get_name(), "PIKAA" );
		TS_ASSERT_EQUALS( commands[2][1]->get_name(), "NATAA" );
		TS_ASSERT_EQUALS( commands[3][1]->get_name(), "PIKAA" );

	}







};
