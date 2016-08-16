// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/devel/cycles/ReindexingMover.cxxtest.hh
/// @brief  
/// @author Kale Kundert

#ifndef INCLUDED_devel_cycles_reindexing_mover_CXXTEST_HH
#define INCLUDED_devel_cycles_reindexing_mover_CXXTEST_HH

#include <cxxtest/TestSuite.h>

#include <devel/cycles/SetupMover.hh>
#include <devel/cycles/ReindexingMover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

using namespace core;
using namespace devel;

class ReindexingMoverTest : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }
	void tearDown() {}

	void test_something() {

		pose::Pose original_pose, reindexed_pose;
		string pdb_path = "devel/cycles/peptides/8.unmarked.pdb";

		import_pose::pose_from_file(original_pose, pdb_path, core::import_pose::PDB_file);
		import_pose::pose_from_file(reindexed_pose, pdb_path, core::import_pose::PDB_file);

		cycles::SetupMover initializer;

		initializer.apply(reindexed_pose);

		original_pose.dump_pdb("original_pose.pdb");
		reindexed_pose.dump_pdb("reindexed_pose.pdb");
	}

};

#endif


