// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/rosettascripts.hh>

// Unit headers
#include <core/pack/task/operation/ClashBasedRepackShell.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// Utility Headers
#include <basic/datacache/DataMap.hh>

using namespace std;
using namespace core::pack::task::operation;
using core::Size;
using core::Real;
using core::import_pose::PDB_file;
using core::import_pose::pose_from_file;
using core::pack::task::TaskFactory;
using core::pack::task::PackerTaskOP;
using core::pose::Pose;

// This class is a pretty thin wrapper around ClashBasedShellSelector, so
// there's no need to extensively test the shell itself.  Really we just need
// to make sure the RosettaScripts interface works.

class ClashBasedRepackShellTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	// Design residue 193 from 2NXX and repack any residues that might interact
	// with it.  I'm intentionally reusing one of the test cases from the
	// ClashBasedShellSelector tests, because we know what it should do.
	void test_xml_tag() {
		Pose pose;
		pose_from_file(pose, "core/pack/task/residue_selector/2NXX.pdb", PDB_file);

		auto allaa_193 = ReadResfileOP( new ReadResfile("core/pack/task/operation/ALLAA_193.resfile") );
		basic::datacache::DataMap data;
		data.add("task_operations", "allaa_193", allaa_193);

		// This is what we're testing.
		auto repack_shell = parse_taskop_tag<ClashBasedRepackShell>(
			R"(<ClashBasedRepackShell task_operations="allaa_193"/>)", data);

		TaskFactory task_factory;
		task_factory.push_back( allaa_193 );
		task_factory.push_back( repack_shell );

		PackerTaskOP task = task_factory.create_task_and_apply_taskoperations( pose );

		// This is the test.  We expect position 193 to design to anything,
		// positions 74,197,423,424,426,427 to repack, and every other position to
		// stay fixed.
		test::UTracer UTR ("core/pack/task/operation/ClashBasedRepackShell.u");
		UTR << *(task);
	}

};
