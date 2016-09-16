// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/cyclic_geometry.cxxtest.hh
/// @brief Unit tests for cyclic peptide pose scoring.
/// @detials Cyclic permutations should score identically.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR("core.scoring.CyclicGeometry_betanov15_Tests.cxxtest");

class CyclicGeometry_betanov15_TwoChainTests : public CxxTest::TestSuite {

public:

	/// @brief Given a pose, remove all disulfides in it.
	///
	void remove_disulfides( core::pose::PoseOP pose ) {
		utility::vector1< std::pair < core::Size, core::Size > > disulfides;
		core::conformation::disulfide_bonds( pose->conformation(), disulfides );

		for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
			core::conformation::break_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second );
		}
	}

	/// @brief Given a pose, add disulfides between the first cysteine and the next.
	///
	void form_disulfides( core::pose::PoseOP pose) {
		bool breaknow(false);
		for ( core::Size ir=1, nres=pose->size(); ir<nres; ++ir ) { //Loop through all residues except the last
			if ( pose->residue(ir).name3() == "CYS" || pose->residue(ir).name3() == "DCS" ) {
				for ( core::Size jr=ir+1; jr<=nres; ++jr ) { //Loop through rest of residues
					if ( pose->residue(jr).name3() == "CYS" || pose->residue(jr).name3() == "DCS" ) {
						core::conformation::form_disulfide( pose->conformation(), ir, jr, true, false );
						breaknow=true;
						break;
					}
				}
			}
			if ( breaknow ) break;
		}
	}

	/// @brief Given a residue type, get its mirror-image type.
	///
	core::chemical::ResidueType const & get_mirror_type( core::chemical::ResidueType const &master_type ) {
		if ( !master_type.is_l_aa() && !master_type.is_d_aa() ) return master_type;
		core::chemical::ResidueTypeSet const & residue_type_set = *master_type.residue_type_set();
		if ( master_type.is_l_aa() ) {
			return residue_type_set.name_map( "D"+master_type.name() );
		}
		return residue_type_set.name_map( master_type.name().substr(1) );
	}

	/// @brief Given a pose, flip the L-residues to D-residues.
	///
	void flip_residues( core::pose::Pose &pose ) {
		// Create the new residue and replace it
		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {
			core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( get_mirror_type( pose.residue(ir).type() ), pose.residue( ir ), pose.conformation());
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( ir ), *new_res, pose.conformation(), true );
			pose.replace_residue( ir, *new_res, false );
		}
	}

	/// @brief Given a pose, construct its mirror image.
	///
	core::pose::PoseOP
	mirror_pose( core::pose::PoseCOP master) {
		core::pose::PoseOP refpose( master->clone() );
		remove_disulfides(refpose);
		core::pose::PoseOP output_pose( refpose->clone() );
		flip_residues(*output_pose);

		for ( core::Size ir=1, irmax=output_pose->size(); ir<=irmax; ++ir ) {
			for ( core::Size ia=1, iamax=output_pose->residue(ir).type().natoms(); ia<=iamax; ++ia ) {
				core::id::AtomID const curatom(ia, ir);
				numeric::xyzVector<core::Real> xyztemp( refpose->xyz( curatom ) );
				xyztemp.z( -1.0*xyztemp.z() );
				output_pose->set_xyz( curatom, xyztemp );
			}
		}

		form_disulfides(output_pose);

		return output_pose;
	}

	/// @brief Given an input cyclic pose, circularly permute by 1 residue.
	core::pose::PoseOP permute(
		core::pose::PoseOP input_pose
	) {
		core::pose::PoseOP ref_pose( input_pose->clone() );
		remove_disulfides(ref_pose);
		core::pose::PoseOP new_pose( new core::pose::Pose );

		for ( core::Size i=2, imax=ref_pose->size(); i<=imax; ++i ) {
			core::conformation::ResidueOP new_rsd( ref_pose->residue(i).clone() );
			if ( i == 2 ) {
				new_pose->append_residue_by_jump(*new_rsd, 1);
			} else {
				new_pose->append_residue_by_bond(*new_rsd, false);
			}
		}
		new_pose->append_residue_by_bond( *(ref_pose->residue(1).clone()), false);
		new_pose->conformation().declare_chemical_bond(1, "N", new_pose->size(), "C");
		form_disulfides(new_pose);
		return new_pose;
	}

	void setUp() {
		core_init_with_additional_options( "-beta_nov15 -score:weights beta_nov15.wts -symmetric_gly_tables true -write_all_connect_info -connect_info_cutoff 0.0" );

		// Pull in the two-chain cyclic peptide pose (12 + 12 = 24 residues)
		core::pose::PoseOP initial_pose_2chain( new core::pose::Pose );
		core::import_pose::pose_from_file( *initial_pose_2chain, "core/scoring/twochain_cycpep.pdb", core::import_pose::PDB_file );

		// Strip off termini and connect the ends:
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::UPPER_TERMINUS_VARIANT, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 12 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::LOWER_TERMINUS_VARIANT, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 13 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::UPPER_TERMINUS_VARIANT, 24 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_LOWER, 24 );
		core::pose::remove_variant_type_from_pose_residue( *initial_pose_2chain, core::chemical::CUTPOINT_UPPER, 24 );
		initial_pose_2chain->conformation().declare_chemical_bond(1, "N", 24, "C");
		initial_pose_2chain->conformation().declare_chemical_bond(13, "N", 12, "C");
		remove_disulfides(initial_pose_2chain);
		form_disulfides(initial_pose_2chain);
		for ( core::Size ir=1, irmax=initial_pose_2chain->size(); ir<=irmax; ++ir ) {
			initial_pose_2chain->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(ir);
		}

		poses_2chain_.push_back(initial_pose_2chain);
		mirror_poses_2chain_.push_back( mirror_pose( poses_2chain_[1] ) );
		for ( core::Size i=1; i<=23; ++i ) {
			poses_2chain_.push_back( permute( poses_2chain_[i] ) );
			mirror_poses_2chain_.push_back( mirror_pose( poses_2chain_[i+1] ) );

			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(1);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(13);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(12);
			poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(24);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(1);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(13);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(12);
			mirror_poses_2chain_[i+1]->conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(24);

		}
		/*initial_pose_2chain->dump_pdb("vtemp00.pdb"); //DELETE ME
		poses_2chain_[1]->dump_pdb("vtemp01.pdb"); //DELETE ME
		poses_2chain_[2]->dump_pdb("vtemp02.pdb"); //DELETE ME
		poses_2chain_[3]->dump_pdb("vtemp03.pdb"); //DELETE ME
		poses_2chain_[4]->dump_pdb("vtemp04.pdb"); //DELETE ME
		poses_2chain_[5]->dump_pdb("vtemp05.pdb"); //DELETE ME
		poses_2chain_[6]->dump_pdb("vtemp06.pdb"); //DELETE ME
		poses_2chain_[7]->dump_pdb("vtemp07.pdb"); //DELETE ME
		poses_2chain_[8]->dump_pdb("vtemp08.pdb"); //DELETE ME
		poses_2chain_[9]->dump_pdb("vtemp09.pdb"); //DELETE ME
		poses_2chain_[10]->dump_pdb("vtemp10.pdb"); //DELETE ME
		poses_2chain_[11]->dump_pdb("vtemp11.pdb"); //DELETE ME
		poses_2chain_[12]->dump_pdb("vtemp12.pdb"); //DELETE ME
		poses_2chain_[13]->dump_pdb("vtemp13.pdb"); //DELETE ME
		poses_2chain_[14]->dump_pdb("vtemp14.pdb"); //DELETE ME
		poses_2chain_[15]->dump_pdb("vtemp15.pdb"); //DELETE ME
		poses_2chain_[16]->dump_pdb("vtemp16.pdb"); //DELETE ME
		poses_2chain_[17]->dump_pdb("vtemp17.pdb"); //DELETE ME
		poses_2chain_[18]->dump_pdb("vtemp18.pdb"); //DELETE ME
		poses_2chain_[19]->dump_pdb("vtemp19.pdb"); //DELETE ME
		poses_2chain_[20]->dump_pdb("vtemp20.pdb"); //DELETE ME
		poses_2chain_[21]->dump_pdb("vtemp21.pdb"); //DELETE ME
		poses_2chain_[22]->dump_pdb("vtemp22.pdb"); //DELETE ME
		poses_2chain_[23]->dump_pdb("vtemp23.pdb"); //DELETE ME
		poses_2chain_[24]->dump_pdb("vtemp24.pdb"); //DELETE ME*/
	}

	void tearDown() {
	}

	/// @brief Test that the same score is returned for all cyclic permutations.
	///
	void cyclic_pose_test( core::scoring::ScoreFunctionOP sfxn ) {
		//Score all of the poses and confirm that they're all equal to the first
		for ( core::Size i=1; i<=24; ++i ) {
			(*sfxn)(*poses_2chain_[i]);
			if ( i>1 ) {
				TR << "\tTesting 2-chain scoring with circular permutation of " << i - 1 << " residues." << std::endl;
				TS_ASSERT_DELTA(poses_2chain_[1]->energies().total_energy(), poses_2chain_[i]->energies().total_energy(), std::abs( std::max(poses_2chain_[1]->energies().total_energy(), poses_2chain_[i]->energies().total_energy())/10000.0 ) );
			}
			//Check mirrored geometry, too:
			TR << "\tTesting 2-chain scoring with circular permutation of " << i - 1 << " residues and mirroring." << std::endl;
			(*sfxn)(*mirror_poses_2chain_[i]);
			TS_ASSERT_DELTA(poses_2chain_[1]->energies().total_energy(), mirror_poses_2chain_[i]->energies().total_energy(), std::abs( std::max(poses_2chain_[1]->energies().total_energy(), mirror_poses_2chain_[i]->energies().total_energy())/10000.0 ) );
			TR.flush();
		}
	}

	/// @brief Tests cyclic permutation scoring with the fa_atr scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_atr() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
		TR << "Testing fa_atr score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
		TR << "Testing fa_rep score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_rep scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_rep() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 1.0 );
		TR << "Testing fa_intra_rep score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_sol scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_sol() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_sol score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_intra_sol_xover4 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_intra_sol_xover4() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_intra_sol_xover4, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing fa_intra_sol_xover4 score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the lk_ball_wtd scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_lk_ball_wtd() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::lk_ball_wtd, 1.0 );
		scorefxn->set_weight( core::scoring::fa_sol, 1.0 );
		TR << "Testing lk_ball_wtd score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_elec scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_elec() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_elec, 1.0 );
		TR << "Testing fa_elec score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the pro_close scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_pro_close() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::pro_close, 1.0 );
		TR << "Testing pro_close score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the dslf_fa13 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_dslf_fa13() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::dslf_fa13, 1.0 );
		TR << "Testing dslf_fa13 score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the hbonds scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_hbonds() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_sc, 1.0 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		TR << "Testing hbonds score terms." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the fa_dun scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_fa_dun() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::fa_dun, 1.0 );
		TR << "Testing fa_dun score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the omega scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_omega() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::omega, 1.0 );
		TR << "Testing omega score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama, 1.0 );
		TR << "Testing rama score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the rama_prepro scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_rama_prepro() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::rama_prepro, 1.0 );
		TR << "Testing rama_prepro score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the p_aa_pp scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_p_aa_pp() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::p_aa_pp, 1.0 );
		TR << "Testing p_aa_pp score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the yhh_planarity scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_yhh_planarity() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->set_weight( core::scoring::yhh_planarity, 1.0 );
		TR << "Testing yhh_planarity score term." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full beta_nov15 scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_beta_nov15() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		scorefxn->add_weights_from_file("beta_nov15.wts");
		TR << "Testing full beta_nov15 score function." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

	/// @brief Tests cyclic permutation scoring with the full scorefunction that's the current default (whatever that is).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_cyclic_permutation_current_default_scorefxn() {
		//Set up the scorefunction
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		TR << "Testing full default score function." << std::endl;
		cyclic_pose_test(scorefxn);
		return;
	}

private:
	utility::vector1 < core::pose::PoseOP > poses_2chain_;
	utility::vector1 < core::pose::PoseOP > mirror_poses_2chain_;


};
