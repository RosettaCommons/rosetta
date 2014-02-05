// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
// Auto-header: duplicate removed #include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

#include <protocols/ligand_docking/RandomConformerMover.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>



static basic::Tracer TR("protocols.ligand_docking.ResidueTorsionRestraintsTest.cxxtest");


class ResidueTorsionRestraintsTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("ZN1")) params_files.push_back("protocols/ligand_docking/ZN1.params");
		if(!residue_set.has_name("CP1")) params_files.push_back("protocols/ligand_docking/7cpa.params");
		residue_set.read_files(params_files,
			ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->element_set( "default" ),
			ChemicalManager::get_instance()->mm_atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->orbital_type_set(FA_STANDARD));//,
			//ChemicalManager::get_instance()->csd_atom_type_set( FA_STANDARD ));
	}

	void tearDown() {}

	void test_randomize_with_constraints() {
		using namespace core::scoring;
		using namespace protocols::ligand_docking;
		using namespace protocols::moves;

		core::Real const score_eps = 1e-6;

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" );

		ScoreFunctionOP sfxn = new ScoreFunction();
		sfxn->set_weight(dihedral_constraint, 1.0);
		TS_ASSERT_DELTA( 0.0, (*sfxn)( pose ), score_eps );

		core::Size const ligres = 309; // the ligand
		TS_ASSERT( pose.residue(ligres).name() == "CP1" );
		TS_ASSERT( !pose.residue(ligres).is_polymer() );

		ResidueTorsionRestraintsOP lig_restraints = new ResidueTorsionRestraints(pose, ligres, 10.0 /*stddev_degrees*/);
		{
			TR << "First trial: randomization with constraints left intact" << std::endl;
			core::Real const start_score = (*sfxn)( pose );
			TR << "Constraint score for original conformation: " << start_score << std::endl;
			TS_ASSERT_DELTA( 0.0, start_score, score_eps );
			MoverOP random_conf = new RandomConformerMover(ligres);
			int score_went_up = 0;
			int const num_trials = 10;
			for(int i = 0; i < num_trials; ++i) {
				// Operating on a copy of the pose is significant because the constraints get cloned.
				core::pose::Pose pose_copy(pose);
				core::Real const old_score = (*sfxn)( pose_copy );
				TS_ASSERT_DELTA( start_score, old_score, score_eps );
				random_conf->apply(pose_copy);
				core::Real const new_score = (*sfxn)( pose_copy );
				TR << "Constraint score for random conformation: " << new_score << std::endl;
				if(new_score > start_score + 1.0) score_went_up += 1;
			}
			TR << "Constraint score increased on " << score_went_up << " out of " << num_trials << " random trials." << std::endl;
			TS_ASSERT( score_went_up > (9*num_trials)/10 );
		}

		{
			TR << "Second trial: randomization with constraints updated" << std::endl;
			core::Real const start_score = (*sfxn)( pose );
			TR << "Constraint score for original conformation: " << start_score << std::endl;
			TS_ASSERT_DELTA( 0.0, start_score, score_eps );
			core::Size const start_num_constraints = pose.constraint_set()->get_all_constraints().size();
			TR << "Pose starts with " << start_num_constraints << " constraints" << std::endl;
			TS_ASSERT( start_num_constraints == 15 );
			MoverOP random_conf = new RandomConformerMover(ligres);
			UnconstrainedTorsionsMover::Restraints restraints;
			restraints.push_back( lig_restraints );
			MoverOP better_random_conf = new UnconstrainedTorsionsMover( random_conf, restraints );
			int score_went_up = 0;
			int const num_trials = 10;
			for(int i = 0; i < num_trials; ++i) {
				// Operating on a copy of the pose is significant because the constraints get cloned.
				core::pose::Pose pose_copy(pose);
				core::Real const old_score = (*sfxn)( pose_copy );
				TS_ASSERT_DELTA( start_score, old_score, score_eps );
				better_random_conf->apply(pose_copy);
				// An old bug in ResidueTorsionRestraints / UnconstrainedTorsionsMover
				// could sometimes lead to doubling up of constraints...
				TR << "Pose now has " << pose_copy.constraint_set()->get_all_constraints().size() << " constraints" << std::endl;
				TS_ASSERT( pose_copy.constraint_set()->get_all_constraints().size() == start_num_constraints );
				core::Real const new_score = (*sfxn)( pose_copy );
				TR << "Constraint score for random conformation: " << new_score << std::endl;
				if(new_score > start_score + score_eps) score_went_up += 1;
			}
			TR << "Constraint score increased on " << score_went_up << " out of " << num_trials << " random trials." << std::endl;
			TS_ASSERT( score_went_up == 0 );
		}
	}
};

