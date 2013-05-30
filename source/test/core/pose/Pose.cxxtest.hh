// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/pose/Pose.cxxtest.hh
/// @brief  unit tests for core::pose::Pose
/// @author

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


namespace test_pose {

// faux observer for Pose observer interface
struct Obs {
	typedef core::Size Size;
	typedef core::pose::signals::DestructionEvent DestructionEvent;
	typedef core::pose::signals::GeneralEvent GeneralEvent;
	typedef core::pose::signals::ConformationEvent ConformationEvent;
	typedef core::pose::signals::EnergyEvent EnergyEvent;

	Obs() : count( 0 ), g_count( 0 ) {}

	void on_destruction_change( DestructionEvent const & ) { ++count; }
	void on_general_change( GeneralEvent const & ) { ++g_count; }
	void on_conformation_change( ConformationEvent const & ) { ++count; }
	void on_energy_change( EnergyEvent const & ) { ++count; }

	int count;
	int g_count;
};

} // namespace test_pose


class PoseTests : public CxxTest::TestSuite
{


public: //setup


	PoseTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // tests


	/// @brief test Pose observer interface
	void test_Pose_observer() {
		using namespace core::scoring;
		using core::pose::Pose;

		using namespace test_pose;

		ScoreFunction scorefxn = *getScoreFunction();

		Pose * pose = new Pose;
		core::import_pose::pose_from_pdb( *pose, "core/pose/pdbinfo_test_in.pdb" );

		// attach observer
		Obs obs;
		pose->attach_destruction_obs( &Obs::on_destruction_change, &obs );
		pose->attach_general_obs( &Obs::on_general_change, &obs );
		pose->attach_conformation_obs( &Obs::on_conformation_change, &obs );
		pose->attach_energy_obs( &Obs::on_energy_change, &obs );

		// test ConformationEvent & GeneralEvent
		pose->set_phi( 3, 45 );
		pose->residue( 2 ); // force refold
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT_EQUALS( obs.g_count, 1 );
		obs.count = 0;
		obs.g_count = 0;

		// test EnergyEvent & GeneralEvent
		scorefxn( *pose );
		TS_ASSERT_EQUALS( obs.count, 1 );
		TS_ASSERT_EQUALS( obs.g_count, 1 );
		obs.count = 0;
		obs.g_count = 0;

		// test DestructionEvent
		delete pose;
		TS_ASSERT_EQUALS( obs.count, 1 );
		obs.count = 0;
	}

	void test_append_pose_by_jump() 
	{
		using namespace core::pose;

		Pose pose;
		make_pose_from_sequence(pose, "TE", "fa_standard");
		pose.pdb_info(PDBInfoOP(new PDBInfo(pose)));

		Pose pose2;
		make_pose_from_sequence(pose2, "ST", "fa_standard");
		pose2.pdb_info(PDBInfoOP(new PDBInfo(pose2)));

		char pose2_chain = 'B';
		pose2.pdb_info()->set_chains(pose2_chain);

		Pose work_pose(pose);
		work_pose.append_pose_by_jump(pose2, 1);

		TS_ASSERT_EQUALS(work_pose.n_residue(), pose.n_residue() + pose2.n_residue());
		TS_ASSERT_EQUALS(work_pose.sequence(), pose.sequence() + pose2.sequence());

		for (core::Size i = 1; i <= pose.n_residue(); i++)
		{
			TS_ASSERT_DELTA(
					work_pose.residue(i).xyz(1),
					pose.residue(i).xyz(1),
					1e-6);
			TS_ASSERT_EQUALS(
					work_pose.pdb_info()->chain(i),
					pose.pdb_info()->chain(i));
		}

		for (core::Size i = 1; i <= pose2.n_residue(); i++)
		{
			TS_ASSERT_DELTA(
					work_pose.residue(i + pose.n_residue()).xyz(1),
					pose2.residue(i).xyz(1),
					1e-6);
			TS_ASSERT_EQUALS(
					work_pose.pdb_info()->chain(i + pose.n_residue()),
					pose2.pdb_info()->chain(i));
		}
	}
};
