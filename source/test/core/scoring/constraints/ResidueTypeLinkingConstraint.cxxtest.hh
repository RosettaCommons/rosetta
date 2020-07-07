// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueTypeLinkingConstraint.cxxtest.hh
/// @brief  test suite for residue type linking constraints
/// @author Alex Sevy
/// @author Brahm Yachnin (brahm.yachnin@visterrainc.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class ResidueTypeLinkingConstraintTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_constraints()
	{
		using namespace core;
		using namespace core::scoring::constraints;

		pose::Pose start_pose( create_test_in_pdb_pose() );
		pose::Pose pose1 (start_pose);
		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 10, 1.0
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 0 );
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 2, 1.0
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 1 );
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 3, 2.0
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 3 );
	}

	// Test ResidueTypeLinkingConstraints if using a std::map
	void test_map_constraints()
	{
		using namespace core;
		using namespace core::chemical;
		using namespace core::scoring::constraints;

		// Set up the pose and scorefxn
		pose::Pose start_pose( create_test_in_pdb_pose() );
		pose::Pose pose1 (start_pose);
		scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared<scoring::ScoreFunction>() );
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);

		// Create a bonus_map to use with the ResidueTypeLinkingConstraint
		// Note: PDB number for this pose is (pose number + 48)
		std::map<std::pair<AA, AA>, core::Real> bonus_map;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("D"), aa_from_one_or_three("I"))] = 1.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("I"), aa_from_one_or_three("D"))] = 2.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("I"), aa_from_one_or_three("I"))] = 3.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("W"), aa_from_one_or_three("Q"))] = 10.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("P"), aa_from_one_or_three("T"))] = -50.0;

		// Add csts between D49 (1) and A50 (2), does not map to any bonus in the bonus_map
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 2,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Add csts between D49 (1) and I51 (3), maps to a bonus of 1.0
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 3,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 1.0 );

		// Add csts between I51 (3) and I53 (5), maps to a bonus of 3.0
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			3, 5,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 4.0 );

		// Add csts between W59 (11) and Q85 (37), maps to a bonus of 10.0
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			11, 37,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 14.0 );

		// Flip order of csts, so that it's between Q85 (37) and W59 (11), so should have no bonus
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			37, 11,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 14.0 );

		// Add csts between P67 (19) and T92 (44), maps to a NEGATIVE bonus of -50.0
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			19, 44,
			bonus_map
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], -36.0 );
	}

	// Test ResidueTypeLinkingConstraints to link residues to specified identities
	void test_linked_identity_constraints()
	{
		using namespace core;
		using namespace core::scoring::constraints;

		// Set up the pose and scorefxn
		pose::Pose start_pose( create_test_in_pdb_pose() );
		pose::Pose pose1 (start_pose);
		pose::Pose pose2 (start_pose);
		scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared<scoring::ScoreFunction>() );
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);

		// Add a constraint based on sequence position and speicified identities (-10.0 if keep current values)
		pose1.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 3,
			core::chemical::aa_from_name("ASP"), core::chemical::aa_from_name("ILE"),
			-10.0
			));
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], -10.0 );

		// Mutate pose1 so that position 3 is a LEU.  -10.0 bonus should go away.
		protocols::simple_moves::MutateResidue mutate_res1( 3, "LEU" );
		mutate_res1.apply(pose1);
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Add a constraint to add a penalty if a particular sequence is observed
		pose2.add_constraint( utility::pointer::make_shared< scoring::constraints::ResidueTypeLinkingConstraint >(
			1, 3,
			core::chemical::aa_from_name("ASP"), core::chemical::aa_from_name("ASP"),
			5.0
			));
		( *scorefxn )( pose2 );
		TS_ASSERT_EQUALS( pose2.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Mutate pose2 so that position 3 is an ASP.  The penalty of 5.0 should be added.
		protocols::simple_moves::MutateResidue mutate_res2( 3, "ASP" );
		mutate_res2.apply(pose2);
		( *scorefxn )( pose2 );
		TS_ASSERT_EQUALS( pose2.energies().total_energies()[ scoring::res_type_linking_constraint ], 5.0 );
	}

	// Test ResidueTypeLinkingConstraints setters
	void test_residue_type_linking_constraints_setters()
	{
		using namespace core;
		using namespace core::chemical;
		using namespace core::scoring::constraints;

		// Set up the pose and scorefxn
		pose::Pose start_pose( create_test_in_pdb_pose() );
		scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared<scoring::ScoreFunction>() );
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);

		// Create an empty bonus_map
		std::map<std::pair<AA, AA>, core::Real> empty_bonus_map;

		// Create a bonus_map to use with the ResidueTypeLinkingConstraint
		// Note: PDB number for this pose is (pose number + 48)
		std::map<std::pair<AA, AA>, core::Real> bonus_map;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("D"), aa_from_one_or_three("I"))] = 1.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("I"), aa_from_one_or_three("D"))] = 2.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("I"), aa_from_one_or_three("I"))] = 3.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("W"), aa_from_one_or_three("Q"))] = 10.0;
		bonus_map[std::pair<AA, AA>(aa_from_one_or_three("P"), aa_from_one_or_three("T"))] = -50.0;

		// Make a pose a score with an empty bonus_map_.  Should score to 0.
		pose::Pose pose1(start_pose);
		ResidueTypeLinkingConstraintOP cst1( utility::pointer::make_shared<ResidueTypeLinkingConstraint>(1, 3, empty_bonus_map ) );
		pose1.add_constraint( cst1 );
		( *scorefxn )( pose1 );
		TS_ASSERT_EQUALS( pose1.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Make a new copy of the pose.  Re-assign the bonus_map to cst1.
		pose::Pose pose1a(start_pose);
		cst1->set_bonus_map(bonus_map);
		pose1a.add_constraint( cst1 );
		( *scorefxn )( pose1a );
		TS_ASSERT_EQUALS( pose1a.energies().total_energies()[ scoring::res_type_linking_constraint ], 1.0 );

		// Create a cst that does not target a current residue pair.
		pose::Pose pose2(start_pose);
		ResidueTypeLinkingConstraintOP cst2(utility::pointer::make_shared<ResidueTypeLinkingConstraint > ( 1, 4, bonus_map ) );
		pose2.add_constraint(cst2);
		( *scorefxn )( pose2 );
		TS_ASSERT_EQUALS( pose2.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Change cst2 to target a current residue pair.
		pose::Pose pose2a(start_pose);
		cst2->set_seqpos2(3);
		pose2a.add_constraint(cst2);
		( *scorefxn )( pose2a );
		TS_ASSERT_EQUALS( pose2a.energies().total_energies()[ scoring::res_type_linking_constraint ], 1.0 );

		// Change cst2 to go off target again.
		pose::Pose pose2b(start_pose);
		cst2->set_seqpos1(2);
		pose2b.add_constraint(cst2);
		( *scorefxn )( pose2b );
		TS_ASSERT_EQUALS( pose2b.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Make a pose to score with cst3.  Start off with a case where residue pair is not targeted.
		pose::Pose pose3(start_pose);
		ResidueTypeLinkingConstraintOP cst3(utility::pointer::make_shared<ResidueTypeLinkingConstraint > ( 1, 4, bonus_map ) );
		pose3.add_constraint(cst3);
		( *scorefxn )( pose3 );
		TS_ASSERT_EQUALS( pose3.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Change cst3 to use a default_bonus_ of -18.
		pose::Pose pose3a(start_pose);
		TS_ASSERT_EQUALS( cst3->get_default_bonus(), 0.0 );
		cst3->set_default_bonus(-18);
		TS_ASSERT_EQUALS( cst3->get_default_bonus(), -18.0 );
		pose3a.add_constraint(cst3);
		( *scorefxn )( pose3a );
		TS_ASSERT_EQUALS( pose3a.energies().total_energies()[ scoring::res_type_linking_constraint ], -18.0 );

		// Change cst3 to go on target.  We won't use default_bonus_ anymore.
		pose::Pose pose3b(start_pose);
		cst3->set_seqpos2(3);
		pose3b.add_constraint(cst3);
		( *scorefxn )( pose3b );
		TS_ASSERT_EQUALS( pose3b.energies().total_energies()[ scoring::res_type_linking_constraint ], 1.0 );
	}

	// Test ResidueTypeLinkingConstraints add_bonus_map_entry functions
	void test_residue_type_linking_constraints_mod_bonus_map()
	{
		using namespace core;
		using namespace core::chemical;
		using namespace core::scoring::constraints;

		// Set up the pose and scorefxn
		pose::Pose start_pose( create_test_in_pdb_pose() );
		scoring::ScoreFunctionOP scorefxn( utility::pointer::make_shared<scoring::ScoreFunction>() );
		scorefxn->set_weight( scoring::res_type_linking_constraint, 1.0);

		// Create a bonus_map to use with the ResidueTypeLinkingConstraint
		// Note: PDB number for this pose is (pose number + 48)
		std::map<std::pair<AA, AA>, core::Real> bonus_map;

		// Create the cst object
		ResidueTypeLinkingConstraintOP cst1( utility::pointer::make_shared<ResidueTypeLinkingConstraint>(1, 3, bonus_map ) );

		// Make two copies of it
		ResidueTypeLinkingConstraintOP cst2( utility::pointer::dynamic_pointer_cast< ResidueTypeLinkingConstraint >(cst1->clone()) );
		ResidueTypeLinkingConstraintOP cst3( utility::pointer::dynamic_pointer_cast< ResidueTypeLinkingConstraint >(cst1->clone()) );

		// Assign cst1 to start_pose and test score
		pose::Pose pose(start_pose);
		pose.add_constraint(cst1);
		( *scorefxn )( pose );
		TS_ASSERT_EQUALS( pose.energies().total_energies()[ scoring::res_type_linking_constraint ], 0.0 );

		// Use the three interfaces, add the same linking cst
		cst1->add_bonus_map_entry(std::pair<AA, AA>(aa_from_oneletter_code('D'), aa_from_oneletter_code('I')), 3.0);
		cst2->add_bonus_map_entry(aa_from_oneletter_code('D'), aa_from_oneletter_code('I'), 3.0);
		cst3->add_bonus_map_entry("DI", 3.0);

		// Rescore pose with the 3.0 key value
		pose = pose::Pose(start_pose);
		pose.add_constraint(cst1);
		( *scorefxn )( pose );
		TS_ASSERT_EQUALS( pose.energies().total_energies()[ scoring::res_type_linking_constraint ], 3.0 );

		// Make sure all bonus maps are the same
		TS_ASSERT_EQUALS(cst1->get_bonus_map(), cst2->get_bonus_map());
		TS_ASSERT_EQUALS(cst1->get_bonus_map(), cst3->get_bonus_map());

		// Check that trying to replace the same entry throws an error
		TS_ASSERT_THROWS(cst2->add_bonus_map_entry(aa_from_oneletter_code('D'), aa_from_oneletter_code('I'), 10.0), utility::excn::Exception & );

		// Check that trying to replace the same entry with overwrite=true works
		TS_ASSERT_THROWS_NOTHING(cst3->add_bonus_map_entry("DI", 10.0, /*overwrite=*/true));

		// Check that trying to add an entry using string syntax that isn't two letters fails
		TS_ASSERT_THROWS(cst2->add_bonus_map_entry("D", 10.0), utility::excn::Exception & );
		TS_ASSERT_THROWS(cst2->add_bonus_map_entry("DIY", 10.0), utility::excn::Exception & );

		// Add another cst1 entry, which should have no effect
		pose = pose::Pose(start_pose);
		cst1->add_bonus_map_entry("DQ", 10.0);
		pose.add_constraint(cst1);
		( *scorefxn )( pose );
		TS_ASSERT_EQUALS( pose.energies().total_energies()[ scoring::res_type_linking_constraint ], 3.0 );

		// Modify the first cst1 entry and make sure score changes
		pose = pose::Pose(start_pose);
		cst1->add_bonus_map_entry("DI", -10.0, /*overwrite=*/true);
		pose.add_constraint(cst1);
		( *scorefxn )( pose );
		TS_ASSERT_EQUALS( pose.energies().total_energies()[ scoring::res_type_linking_constraint ], -10.0 );
	}

	void test_serialize_ResidueTypeLinkingConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		core::pose::Pose pose1( create_test_in_pdb_pose() );
		ConstraintOP instance( utility::pointer::make_shared<ResidueTypeLinkingConstraint>( 1, 10, 1.0 ) ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a ResidueTypeLinkingConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< ResidueTypeLinkingConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
