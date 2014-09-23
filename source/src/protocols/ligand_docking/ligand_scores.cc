// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/MultiResidueLigandDock.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/ligand_scores.hh>
#include <core/pose/util.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>

//#include <protocols/qsar/qsarTypeManager.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>

#include <protocols/jd2/Job.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <ObjexxFCL/FArray1.io.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/MoveMap.fwd.hh> // needed?

#include <sstream>

#include <utility/vector1.hh>

//Auto Headers
namespace protocols {
namespace ligand_docking {

using namespace ObjexxFCL;

///@brief append interface_delta scores
void
append_interface_deltas(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		const core::scoring::ScoreFunctionOP scorefxn,
		std::string const & prefix
){
	core::pose::PoseOP after_copy( new core::pose::Pose( after ) );

	char const ligand_chain= core::pose::get_chain_from_jump_id(jump_id, after);

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	core::Real const together_score = (*scorefxn)( *after_copy );
	core::scoring::EnergyMap const together_energies = after_copy->energies().total_energies();
	core::Real const initial_fa_rep = after_copy->energies().total_energies()[ core::scoring::fa_rep ];
	protocols::rigid::RigidBodyTransMover trans_mover( *after_copy, jump_id );
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move together
	trans_mover.step_size(1);
	trans_mover.apply( *after_copy );
	(*scorefxn)( *after_copy );
	core::Real const push_together_fa_rep = after_copy->energies().total_energies()[ core::scoring::fa_rep ];
	bool const are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);

	{
		std::ostringstream touching;
		touching << "ligand_is_touching_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(touching.str(), are_touching);
		}else
		{
			job->add_string_real_pair(prefix + "_" + touching.str(), are_touching);
		}
	}

	// Now pull apart by 500 A to determine the reference E for calculating interface E.
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move apart
	trans_mover.step_size(500); // make sure they're fully separated!
	trans_mover.apply( *after_copy );
	core::Real const separated_score = (*scorefxn)( *after_copy );
	core::scoring::EnergyMap const separated_energies = after_copy->energies().total_energies();

	{
		std::ostringstream delta;
		delta<< "interface_delta_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(delta.str(), together_score - separated_score);
		}else
		{
			job->add_string_real_pair(prefix + "_" + delta.str(), together_score - separated_score);
		}
	}

	// Interface delta, broken down by component
	for(int i = 1; i <= core::scoring::n_score_types; ++i) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);

		if ( !scorefxn->has_nonzero_weight(ii) ) continue;
		core::Real component_score= scorefxn->get_weight(ii) * (together_energies[ii] - separated_energies[ii]);
		{
			std::ostringstream if_score;
			if_score<< "if_"<< ligand_chain<< '_' << name_from_score_type(ii);
			if(prefix == "")
			{
				job->add_string_real_pair(if_score.str(), component_score);
			}else
			{
				job->add_string_real_pair(prefix + "_" + if_score.str(), component_score);
			}

		}
	}
}

void
append_interface_deltas(
	core::Size jump_id,
	protocols::jd2::JobOP job,
	core::pose::Pose const & after,
	const core::scoring::ScoreFunctionOP scorefxn,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
	)
{
	core::pose::PoseOP after_copy( new core::pose::Pose( after ) );

	char const ligand_chain= core::pose::get_chain_from_jump_id(jump_id, after);
	core::Size chain_id = core::pose::get_chain_id_from_jump_id(jump_id, after);
	core::conformation::ResidueCOPs residues(core::pose::get_chain_residues(after, chain_id));

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	core::Real const together_score = (*normalization_function)((*scorefxn)( *after_copy ),residues);
	core::scoring::EnergyMap const together_energies = after_copy->energies().total_energies();
	core::Real const initial_fa_rep = (*normalization_function)(after_copy->energies().total_energies()[ core::scoring::fa_rep ],residues);
	protocols::rigid::RigidBodyTransMover trans_mover( *after_copy, jump_id );
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move together
	trans_mover.step_size(1);
	trans_mover.apply( *after_copy );
	(*scorefxn)( *after_copy );
	core::Real const push_together_fa_rep = (*normalization_function)(after_copy->energies().total_energies()[ core::scoring::fa_rep ],residues);
	bool const are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);

	{
		std::ostringstream touching;
		touching << "ligand_is_touching_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(touching.str(), are_touching);
		}else
		{
			job->add_string_real_pair(prefix + "_" + touching.str(), are_touching);
		}
	}

	// Now pull apart by 500 A to determine the reference E for calculating interface E.
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move apart
	trans_mover.step_size(500); // make sure they're fully separated!
	trans_mover.apply( *after_copy );
	core::Real const separated_score = (*normalization_function)((*scorefxn)( *after_copy ),residues);
	core::scoring::EnergyMap const separated_energies = after_copy->energies().total_energies();

	{
		std::ostringstream delta;
		delta<< "interface_delta_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(delta.str(), together_score - separated_score);
		}else
		{
			job->add_string_real_pair(prefix + "_" + delta.str(), together_score - separated_score);
		}

	}

	// Interface delta, broken down by component
	for(int i = 1; i <= core::scoring::n_score_types; ++i) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);

		if ( !scorefxn->has_nonzero_weight(ii) ) continue;

		core::Real component_score= (*normalization_function)(scorefxn->get_weight(ii) * (together_energies[ii] - separated_energies[ii]),residues);
		{
			std::ostringstream if_score;
			if_score<< "if_"<< ligand_chain<< '_' << name_from_score_type(ii);
			if(prefix == "")
			{
				job->add_string_real_pair(if_score.str(), component_score);
			}else
			{
				job->add_string_real_pair(prefix + "_" + if_score.str(), component_score);
			}
		}
	}
}

///@brief Another interesting metric -- how far does the ligand centroid move?
///@brief Large values indicate we're outside of the intended binding site.
void
append_ligand_travel(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
){
	core::Vector const before_vector= protocols::geometry::downstream_centroid_by_jump(before, jump_id);
	core::Vector const after_vector= protocols::geometry::downstream_centroid_by_jump(after, jump_id);

	char const ligand_chain= core::pose::get_chain_from_jump_id(jump_id, after);
	{
		std::ostringstream centroid;
		centroid<< "ligand_centroid_travel_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(centroid.str(), before_vector.distance(after_vector));
		}else
		{
			job->add_string_real_pair(prefix + "_" + centroid.str(), before_vector.distance(after_vector));
		}
	}
}

void append_ligand_grid_scores(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		std::string const & prefix
)
{
	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();

	if (grid_manager->size()==0){
		ligand_scores_tracer << "skipping 'append ligand grid scores'. No grids used." << std::endl;
		return;
	}

	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(after,jump_id));
	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(after,center);
	core::Size const chain_id = core::pose::get_chain_id_from_jump_id(jump_id,after);
	core::Real total_score = grid_manager->total_score(after,chain_id);
	char const ligand_chain=core::pose::get_chain_from_jump_id(jump_id,after);

	qsar::scoring_grid::ScoreMap grid_scores(grid_manager->get_cached_scores());
	BOOST_FOREACH(qsar::scoring_grid::ScoreMap::value_type grid_score, grid_scores){
		std::ostringstream score_label;
		score_label << grid_score.first << "_grid_" <<ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(score_label.str(),grid_score.second);
		}else
		{
			job->add_string_real_pair(prefix + "_" + score_label.str(),grid_score.second);
		}

	}

	std::ostringstream score_label;
	score_label << "total_score_" <<ligand_chain;
	if(prefix == "")
	{
		job->add_string_real_pair(score_label.str(),total_score);
	}else
	{
		job->add_string_real_pair(prefix + "_" + score_label.str(),total_score);
	}
}

void append_ligand_grid_scores(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		std::string const & prefix,
		protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
)
{
	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();

	if (grid_manager->size()==0){
		ligand_scores_tracer << "skipping 'append ligand grid scores'. No grids used.";
		return;
	}

	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(after,jump_id));
	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(after,center);
	core::Size const chain_id = core::pose::get_chain_id_from_jump_id(jump_id,after);
	core::conformation::ResidueCOPs residues(core::pose::get_chain_residues(after, chain_id));
	core::Real total_score = (*normalization_function)(grid_manager->total_score(after,chain_id),residues);
	char const ligand_chain=core::pose::get_chain_from_jump_id(jump_id,after);

	qsar::scoring_grid::ScoreMap grid_scores(grid_manager->get_cached_scores());
	BOOST_FOREACH(qsar::scoring_grid::ScoreMap::value_type grid_score, grid_scores){
		std::ostringstream score_label;
		score_label << grid_score.first << "_grid_" <<ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(score_label.str(),(*normalization_function)(grid_score.second,residues));
		}else
		{
			job->add_string_real_pair(prefix + "_" + score_label.str(),(*normalization_function)(grid_score.second,residues));
		}
	}

	std::ostringstream score_label;
	score_label << "total_score_" <<ligand_chain;
	if(prefix == "")
	{
		job->add_string_real_pair(score_label.str(),total_score);
	}else
	{
		job->add_string_real_pair(prefix + "_" + score_label.str(),total_score);
	}

}

///@brief Calculate radius of gyration for downstream non-H atoms
///@brief Ligands tend to bind in outstretched conformations...
void
append_radius_of_gyration(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		std::string const & prefix
){

	core::Real lig_rg = 0;
	int lig_rg_natoms = 0;
	core::Vector const before_vector= protocols::geometry::downstream_centroid_by_jump(before, jump_id);
	FArray1D_bool is_upstream ( before.total_residue(), false );
	before.fold_tree().partition_by_jump( jump_id, is_upstream );
	for(core::Size i = 1, i_end = before.total_residue(); i <= i_end; ++i) {
		if( is_upstream(i) ) continue; // only downstream residues
		core::conformation::Residue const & rsd = before.residue(i);
		for(core::Size j = 1, j_end = rsd.nheavyatoms(); j <= j_end; ++j) {
			lig_rg += before_vector.distance_squared( rsd.xyz(j) );
			lig_rg_natoms += 1;
		}
	}
	lig_rg = std::sqrt( lig_rg / lig_rg_natoms );

	char const ligand_chain= core::pose::get_chain_from_jump_id(jump_id, before);
	{
		std::ostringstream centroid;
		centroid<< "ligand_radius_of_gyration_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(centroid.str(), lig_rg);
		}else
		{
			job->add_string_real_pair(prefix + "_" + centroid.str(), lig_rg);
		}
	}


}

void
append_ligand_RMSD(
		core::Size const jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
){
	assert(before.num_jump() >= jump_id);
	assert(after.num_jump() >= jump_id);

	core::Size chain_id= core::pose::get_chain_id_from_jump_id(jump_id, before);
	core::Size const begin = before.conformation().chain_begin(chain_id);
	core::Size const end = before.conformation().chain_end(chain_id);

	if (end-begin > 0){
		append_multi_residue_ligand_RMSD(jump_id, job, before, after,prefix);
	}else{
		append_automorphic_rmsd(begin, job, before, after,prefix);
	}
}

void
append_multi_residue_ligand_RMSD(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
){
	core::pose::Pose before_ligand;
	core::pose::Pose after_ligand;
	{ ///TODO make this section a function
		core::Size const before_first_residue = before.fold_tree().downstream_jump_residue(jump_id);
		core::Size const after_first_residue = after.fold_tree().downstream_jump_residue(jump_id);
		core::Size const before_chain= before.chain(before_first_residue);
		core::Size const after_chain= before.chain(after_first_residue);
		core::pose::Pose before_copy(before);
		core::pose::PoseOP after_copy( new core::pose::Pose(after) );
		before_copy.remove_constraints();/// TODO fix split_by_chain to avoid this
		after_copy->remove_constraints();/// TODO fix split_by_chain to avoid this
		before_ligand= *before_copy.split_by_chain(before_chain);
		after_ligand= *after_copy->split_by_chain(after_chain);
	}

	char const ligand_chain= core::pose::get_chain_from_jump_id(jump_id, after);
	{
		core::Real ligand_rms_no_super= core::scoring::rmsd_no_super(
				before_ligand,
				after_ligand,
				core::scoring::is_ligand_heavyatom
		);
		std::ostringstream rms_no_super;
		rms_no_super<< "ligand_rms_no_super_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(rms_no_super.str(), ligand_rms_no_super);
		}else
		{
			job->add_string_real_pair(prefix + "_" + rms_no_super.str(), ligand_rms_no_super);
		}
	}
	{
		core::Real ligand_rms_with_super= core::scoring::rmsd_with_super(
				before_ligand,
				after_ligand,
				core::scoring::is_ligand_heavyatom
		);
		std::ostringstream rms_with_super;
		rms_with_super<< "ligand_rms_with_super_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(rms_with_super.str(), ligand_rms_with_super);
		}else
		{
			job->add_string_real_pair(prefix + "_" + rms_with_super.str(), ligand_rms_with_super);
		}
	}
}

void
append_automorphic_rmsd(
		core::Size ligand_residue_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
){

	char const ligand_chain= before.pdb_info()->chain(ligand_residue_id);
	{
		core::Real ligand_rms_no_super= core::scoring::automorphic_rmsd(
				before.residue(ligand_residue_id),
				after.residue(ligand_residue_id),
				false /*don't superimpose*/
		);

		std::ostringstream rms_no_super;
		rms_no_super<< "ligand_rms_no_super_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(rms_no_super.str(), ligand_rms_no_super);
		}else
		{
			job->add_string_real_pair(prefix + "_" + rms_no_super.str(), ligand_rms_no_super);
		}
	}
	{
		core::Real ligand_rms_with_super= core::scoring::automorphic_rmsd(
				before.residue(ligand_residue_id),
				after.residue(ligand_residue_id),
				true /*superimpose*/
		);
		std::ostringstream rms_with_super;
		rms_with_super<< "ligand_rms_with_super_"<< ligand_chain;
		if(prefix == "")
		{
			job->add_string_real_pair(rms_with_super.str(), ligand_rms_with_super);
		}else
		{
			job->add_string_real_pair(prefix + "_" + rms_with_super.str(), ligand_rms_with_super);
		}
	}

}

} // namespace ligand_docking
} // namespace protocols
