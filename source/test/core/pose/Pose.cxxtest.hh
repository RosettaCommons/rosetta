// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
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
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


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

		ScoreFunctionOP scorefxn = get_score_function();

		Pose * pose = new Pose;
		core::import_pose::pose_from_file( *pose, "core/pose/pdbinfo_test_in.pdb" , core::import_pose::PDB_file);

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
		scorefxn->score( *pose );
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
		pose.pdb_info(PDBInfoOP( new PDBInfo(pose) ));

		Pose pose2;
		make_pose_from_sequence(pose2, "ST", "fa_standard");
		pose2.pdb_info(PDBInfoOP( new PDBInfo(pose2) ));

		char pose2_chain = 'B';
		pose2.pdb_info()->set_chains(pose2_chain);

		Pose work_pose(pose);
		work_pose.append_pose_by_jump(pose2, 1);

		TS_ASSERT_EQUALS(work_pose.n_residue(), pose.n_residue() + pose2.n_residue());
		TS_ASSERT_EQUALS(work_pose.sequence(), pose.sequence() + pose2.sequence());

		for ( core::Size i = 1; i <= pose.n_residue(); i++ ) {
			TS_ASSERT_DELTA(
				work_pose.residue(i).xyz(1),
				pose.residue(i).xyz(1),
				1e-6);
			TS_ASSERT_EQUALS(
				work_pose.pdb_info()->chain(i),
				pose.pdb_info()->chain(i));
		}

		for ( core::Size i = 1; i <= pose2.n_residue(); i++ ) {
			TS_ASSERT_DELTA(
				work_pose.residue(i + pose.n_residue()).xyz(1),
				pose2.residue(i).xyz(1),
				1e-6);
			TS_ASSERT_EQUALS(
				work_pose.pdb_info()->chain(i + pose.n_residue()),
				pose2.pdb_info()->chain(i));
		}
	}

	void test_update_residue_neighbors() {
		// These are tests mainly just to make sure that the calls don't crash.
		using namespace core::pose;

		Pose empty;
		empty.update_residue_neighbors();

		Pose pose;
		core::import_pose::pose_from_file( pose, "core/pose/pdbinfo_test_in.pdb" , core::import_pose::PDB_file);
		pose.update_residue_neighbors();
	}

	void test_serialize_pose() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using namespace core::conformation;
		using namespace core::pose;
		using namespace core::id;
		using namespace core::kinematics;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;

		// TO DO: Test constraint set propertly observes the deserialized pose (i.e. following residue insertion)
		// TO DO: Test CacheableObservers properly observe the deserialized pose
		// TO DO: Test that a symmetric conformation is properly deserialized inside a pose

		PoseOP trpcage = create_trpcage_ideal_poseop();
		HarmonicFuncOP harm_0_1( new HarmonicFunc( 0, 1 ) );
		ConstraintOP apc_1_10( new AtomPairConstraint( AtomID(3,1), AtomID(3,10), harm_0_1 ));
		ConstraintOP apc_1_11( new AtomPairConstraint( AtomID(3,1), AtomID(3,11), harm_0_1 ));
		trpcage->add_constraint( apc_1_10 );

		// Now serialize the Pose
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( trpcage );
		}

		PoseOP trpcage_copy;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( trpcage_copy );
		}

		TS_ASSERT( trpcage->total_residue() == trpcage_copy->total_residue() );

		// Make sure constraint set was serialized properly
		TS_ASSERT( ! trpcage_copy->constraint_set()->is_empty() );
		TS_ASSERT( trpcage_copy->constraint_set()->get_all_constraints().size() == 1 );
		TS_ASSERT( *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_10 );
		// make sure that constraint set is observing the conformation
		// do something silly by making another copy of residue 15 and inserting it between
		// residue 5 and residue 6 (set a foldtree that allows this insertion first)
		FoldTree ft;
		ft.add_edge( 1, 5, Edge::PEPTIDE );  ft.add_edge( 5, 10, 1  );
		ft.add_edge( 10, 6, Edge::PEPTIDE ); ft.add_edge( 10, 20, Edge::PEPTIDE );
		trpcage_copy->fold_tree( ft );
		trpcage_copy->insert_residue_by_bond( trpcage_copy->residue(15), 6, 5 );
		TS_ASSERT( ! ( *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_10 ));
		TS_ASSERT(     *trpcage_copy->constraint_set()->get_all_constraints()[1] == *apc_1_11 );

		// Make sure that the deserialized PDBInfo object is reestablished as observing the
		// Pose's Conformation
		core::conformation::ConformationCOP obs_by_pdb_info = trpcage_copy->pdb_info()->is_observing().lock();
		TS_ASSERT_EQUALS( obs_by_pdb_info.get(), &trpcage_copy->conformation() );


#endif
	}

};
