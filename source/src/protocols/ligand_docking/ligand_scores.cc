// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief
/// @author Gordon Lemmon
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/ligand_docking/ligand_scores.hh>
#include <core/pose/util.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>

#include <core/select/residue_selector/ChainSelector.hh>
#include <core/kinematics/util.hh>

#include <core/types.hh>
#include <core/scoring/ScoreType.hh>

#include <basic/Tracer.hh>

#include <core/kinematics/MoveMap.fwd.hh> // needed?

#include <sstream>

#include <utility/vector1.hh>

//Auto Headers
namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer TR("protocols.ligand_docking.ligand_scores");

using namespace ObjexxFCL;

std::map< std::string, core::Real >
get_interface_deltas(
	char chain,
	core::pose::Pose const & after,
	const core::scoring::ScoreFunctionOP scorefxn,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
) {
	std::map< std::string, core::Real > retval;

	core::pose::PoseOP after_copy( new core::pose::Pose( after ) );

	core::select::residue_selector::ChainSelector chain_select( chain );
	utility::vector1< bool > chain_residues( chain_select.apply( *after_copy ) );

	core::Size jump_id( core::kinematics::jump_which_partitions( after_copy->fold_tree(), chain_residues ) );

	if ( jump_id == 0 ) {
		TR.Warning << "In get_interface_deltas(), passed pose is not set up for proper docking of chain " << chain << " reorganizing to make it so." << std::endl;
		core::kinematics::FoldTree new_fold_tree( core::kinematics::get_foldtree_which_partitions( after_copy->fold_tree(), chain_residues ) );
		after_copy->fold_tree( new_fold_tree );
		jump_id = core::kinematics::jump_which_partitions( after_copy->fold_tree(), chain_residues );
		debug_assert( jump_id != 0 );
	}

	core::conformation::ResidueCOPs residues;
	if ( normalization_function ) {
		for ( core::Size ii(1); ii <= chain_residues.size(); ++ii ) {
			residues.push_back( after_copy->residue( ii ).get_self_ptr() );
		}
	}

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	core::Real together_score = (*scorefxn)( *after_copy );
	core::scoring::EnergyMap together_energies = after_copy->energies().total_energies(); // Make a copy
	core::Real initial_fa_rep = after_copy->energies().total_energies()[ core::scoring::fa_rep ];
	if ( normalization_function ) {
		together_score = (*normalization_function)( together_score, residues );
		initial_fa_rep = (*normalization_function)( initial_fa_rep, residues );
	}
	protocols::rigid::RigidBodyTransMover trans_mover( *after_copy, jump_id );
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move together
	trans_mover.step_size(1);
	trans_mover.apply( *after_copy );
	(*scorefxn)( *after_copy );
	core::Real push_together_fa_rep = after_copy->energies().total_energies()[ core::scoring::fa_rep ];
	if ( normalization_function ) {
		push_together_fa_rep = (*normalization_function)( push_together_fa_rep, residues );
	}
	bool const are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);

	std::string touching_label( "ligand_is_touching_" );
	touching_label += chain;
	if ( prefix != "" ) {
		touching_label = prefix + "_" + touching_label;
	}
	retval[ touching_label ] = are_touching;

	// Now pull apart by 500 A to determine the reference E for calculating interface E.
	trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move apart
	trans_mover.step_size(500); // make sure they're fully separated!
	trans_mover.apply( *after_copy );
	core::Real separated_score = (*scorefxn)( *after_copy );
	if ( normalization_function ) {
		separated_score = (*normalization_function)( separated_score, residues );
	}
	core::scoring::EnergyMap const & separated_energies = after_copy->energies().total_energies(); // reference is fine

	std::string delta_label( "interface_delta_" );
	delta_label += chain;
	if ( prefix != "" ) {
		delta_label = prefix + "_" + delta_label;
	}
	retval[ delta_label ] = together_score - separated_score;

	// Interface delta, broken down by component
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);

		if ( !scorefxn->has_nonzero_weight(ii) ) continue;

		core::Real component_score = scorefxn->get_weight(ii) * (together_energies[ii] - separated_energies[ii]);
		if ( normalization_function ) {
			component_score = (*normalization_function)( component_score, residues );
		}
		std::string component_label( "if_" );
		component_label += chain;
		component_label += '_';
		component_label += name_from_score_type(ii);
		if ( prefix != "" ) {
			component_label = prefix + "_" + component_label;
		}
		retval[ component_label ] = component_score;
	}

	return retval;
}

std::map< std::string, core::Real >
get_ligand_travel(
	char chain,
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	std::string const & prefix,
	bool use_ensemble_best /*= false */
) {

	utility::vector1< core::Size > test_resi( core::pose::get_resnums_for_chain( test_pose, chain ) );
	utility::vector1< core::Size > ref_resi( core::pose::get_resnums_for_chain( ref_pose, chain ) );

	if ( test_resi.empty() ) {
		utility_exit_with_message( "In get_ligand_travel: test pose does not have a chain " + utility::to_string( chain ) );
	}
	if ( ref_resi.empty() ) {
		utility_exit_with_message( "In get_ligand_travel: reference pose does not have a chain " + utility::to_string( chain ) );
	}

	core::Real ligand_travel = 0;
	if ( use_ensemble_best ) {
		if ( test_resi.size() != 1 ) {
			utility_exit_with_message( "In get_ligand_travel: cannot use ensemble best method with multi-residue working chain." );
		}
		ligand_travel = get_ligand_travel_ensemble_best( test_pose, ref_pose, test_resi[1], ref_resi );
	} else {
		core::Vector test_center( core::pose::all_atom_center( test_pose, test_resi ) );
		core::Vector ref_center( core::pose::all_atom_center( ref_pose, ref_resi ) );
		ligand_travel = ref_center.distance( test_center );
	}

	std::string label("ligand_centroid_travel_");
	label += chain;
	if ( prefix != "" ) {
		label = prefix + "_" + label;
	}

	std::map< std::string, core::Real > retval;
	retval[ label ] = ligand_travel;
	return retval;
}

core::Real
get_ligand_travel_ensemble_best(
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	core::Size test_residue_id,
	utility::vector1< core::Size > const & ref_residue_ids
) {

	core::Real min_travel = 999999999;
	core::Vector test_center( core::conformation::all_atom_center( test_pose.residue( test_residue_id ) ) );
	for ( core::Size const & ref_residue_id: ref_residue_ids ) {
		core::Vector ref_center( core::conformation::all_atom_center( ref_pose.residue( ref_residue_id ) ) );
		core::Real travel( ref_center.distance( test_center ) );
		if ( travel < min_travel ) {
			min_travel = travel;
		}
	}

	return min_travel;
}

std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	char chain,
	core::pose::Pose const & test_pose,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
) {
	utility::vector1< core::Size > test_resi( core::pose::get_resnums_for_chain( test_pose, chain ) );

	if ( test_resi.empty() ) {
		utility_exit_with_message( "In get_ligand_grid_scores: test pose does not have a chain " + utility::to_string( chain ) );
	}

	return get_ligand_grid_scores( grid_set_prototype, test_resi, utility::to_string( chain ), test_pose, prefix, normalization_function );
}

std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	core::Size jump_id,
	core::pose::Pose const & test_pose,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
) {

	if ( jump_id == 0 || jump_id > test_pose.num_jump() ) {
		utility_exit_with_message( "In get_ligand_grid_scores: test pose does not have jump_id " + utility::to_string( jump_id ) );
	}

	utility::vector1< core::Size > test_resi( core::kinematics::residues_downstream_of_jump( test_pose.fold_tree(), jump_id ) );

	char const chain =core::pose::get_chain_from_jump_id( jump_id, test_pose );

	return get_ligand_grid_scores( grid_set_prototype, test_resi, utility::to_string( chain ), test_pose, prefix, normalization_function );
}

std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	utility::vector1< core::Size > const & test_resi,
	std::string const & chain,
	core::pose::Pose const & test_pose,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
) {

	std::map< std::string, core::Real > retval;

	core::Vector center( core::pose::all_atom_center( test_pose, test_resi ) );

	if ( grid_set_prototype.size()==0 ) {
		TR << "skipping 'append ligand grid scores'. No grids used." << std::endl;
		return retval;
	}

	qsar::scoring_grid::GridSetCOP grid_set = qsar::scoring_grid::GridManager::get_instance()->get_grids( grid_set_prototype, test_pose, center, chain );

	if ( grid_set->is_normalization_enabled() ) {
		// Note that the GridSet::total_score() function would also normalize, if the normalization function was set in the GridSet
		normalization_function = nullptr; // Not passed by reference, so wouldn't change any upstream pointers
	}

	core::conformation::ResidueCOPs residues;
	if ( normalization_function ) {
		for ( core::Size resi: test_resi ) {
			residues.push_back( test_pose.residue( resi ).get_self_ptr() );
		}
	}

	core::Real total_score( grid_set->total_score( test_pose, test_resi ) );
	if ( normalization_function ) {
		total_score = (*normalization_function)( total_score, residues );
	}

	qsar::scoring_grid::GridSet::ScoreMap grid_scores( grid_set->grid_scores( test_pose, test_resi ) );
	for ( auto const & grid_score : grid_scores ) {
		core::Real component_score( grid_score.second );
		if (  normalization_function ) {
			component_score = (*normalization_function)( component_score, residues );
		}

		std::string label(grid_score.first);
		label += "_grid_" + chain;
		if ( prefix != "" ) {
			label = prefix + "_" + label;
		}

		retval[ label ] = component_score;
	}

	std::string total_label("total_score_" + chain);
	if ( prefix != "" ) {
		total_label = prefix + "_" + total_label;
	}

	retval[ total_label ] = total_score;

	return retval;
}


/// @brief Calculate radius of gyration for downstream non-H atoms
/// @brief Ligands tend to bind in outstretched conformations...
std::map< std::string, core::Real >
get_radius_of_gyration(
	char chain,
	core::pose::Pose const & test_pose,
	std::string const & prefix
) {

	utility::vector1< core::Size > test_resi( core::pose::get_resnums_for_chain( test_pose, chain ) );

	if ( test_resi.empty() ) {
		utility_exit_with_message( "In get_radius_of_gyration: test pose does not have a chain " + utility::to_string( chain ) );
	}

	core::Vector test_center( core::pose::all_atom_center( test_pose, test_resi ) );

	core::Real lig_rg = 0;
	int lig_rg_natoms = 0;
	for ( core::Size resi: test_resi ) {
		core::conformation::Residue const & rsd( test_pose.residue(resi) );
		for ( core::Size jj = 1, jj_end = rsd.nheavyatoms(); jj <= jj_end; ++jj ) {
			lig_rg += test_center.distance_squared( rsd.xyz(jj) );
			lig_rg_natoms += 1;
		}
	}
	lig_rg = std::sqrt( lig_rg / lig_rg_natoms );

	std::string label("ligand_radius_of_gyration_");
	label += chain;
	if ( prefix != "" ) {
		label = prefix + "_" + label;
	}

	std::map< std::string, core::Real > retval;
	retval[ label ] = lig_rg;
	return retval;
}


std::map< std::string, core::Real >
get_ligand_RMSDs(
	char chain,
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	std::string const & prefix,
	bool use_ensemble_best
){

	utility::vector1<core::Size> ref_resi( core::pose::get_resnums_for_chain( ref_pose, chain ) );
	utility::vector1<core::Size> test_resi( core::pose::get_resnums_for_chain( test_pose, chain ) );

	if ( test_resi.empty() ) {
		utility_exit_with_message( std::string("get_ligand_RMSD: Cannot find chain ") + chain + " in test pose.");
	}
	if ( ref_resi.empty() ) {
		utility_exit_with_message( std::string("get_ligand_RMSD: Cannot find chain ") + chain + " in reference pose.");
	}

	if ( use_ensemble_best ) {
		if ( test_resi.size() != 1 ) {
			utility_exit_with_message( "get_ligand_RMSD: Cannot use ensemble best method with multi-residue working chain." );
		}
		return get_automorphic_RMSDs( test_pose, ref_pose, test_resi[1], ref_resi, prefix );
	}

	if ( test_resi.size() != ref_resi.size() ) {
		TR.Error << "In get_ligand_RMSD: test and reference poses have a different number of residues for chain " << chain << std::endl;
		TR.Error << "  test pose has " << test_resi.size() << " residues to the reference pose's " << ref_resi.size() << std::endl;
		utility_exit_with_message("Cannot compute (non-ensemble-best) ligand RMSD - different number of residues in the two structures.");
	}

	if ( test_resi.size() == 1 && ref_resi.size() == 1 ) {
		return get_automorphic_RMSDs( test_pose, ref_pose, test_resi[1], ref_resi, prefix );
	}

	TR.Debug << "get_ligand_RMSDs: Chain " << chain << " has multiple (" << test_resi.size() << ") residues in each pose. Using multi-residue version of RMSD calculation." << std::endl;
	return get_multi_residue_ligand_RMSDs( test_pose, ref_pose, test_resi, ref_resi, chain, prefix );
}

std::map< std::string, core::Real >
get_automorphic_RMSDs(
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	core::Size test_residue_id,
	utility::vector1<core::Size> const & ref_residue_ids,
	std::string const & prefix
){
	debug_assert( ! ref_residue_ids.empty() );

	std::map< std::string, core::Real > retval;

	char const ligand_chain = test_pose.pdb_info()->chain(test_residue_id);
	{
		core::Real ligand_rms_no_super = 9999999999;
		for ( core::Size ref_residue_id: ref_residue_ids ) {
			core::Real current_rmsd = core::scoring::automorphic_rmsd(
				ref_pose.residue(ref_residue_id),
				test_pose.residue(test_residue_id),
				false /*don't superimpose*/
			);
			if ( current_rmsd < ligand_rms_no_super ) {
				ligand_rms_no_super = current_rmsd;
			}
		}

		std::string no_super_label("ligand_rms_no_super_");
		no_super_label += ligand_chain;
		if ( prefix != "" ) {
			no_super_label = prefix + "_" + no_super_label;
		}

		retval[ no_super_label ] = ligand_rms_no_super;
	}
	{
		core::Real ligand_rms_with_super = 9999999999;
		for ( core::Size ref_residue_id: ref_residue_ids ) {
			core::Real current_rmsd = core::scoring::automorphic_rmsd(
				ref_pose.residue(ref_residue_id),
				test_pose.residue(test_residue_id),
				true /*superimpose*/
			);
			if ( current_rmsd < ligand_rms_with_super ) {
				ligand_rms_with_super = current_rmsd;
			}
		}

		std::string with_super_label("ligand_rms_with_super_");
		with_super_label += ligand_chain;
		if ( prefix != "" ) {
			with_super_label = prefix + "_" + with_super_label;
		}

		retval[ with_super_label ] = ligand_rms_with_super;
	}

	return retval;
}

std::map< std::string, core::Real >
get_multi_residue_ligand_RMSDs(
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & test_residue_ids,
	utility::vector1< core::Size > const & ref_residue_ids,
	char chain,
	std::string const & prefix
) {
	debug_assert( test_residue_ids.size() == ref_residue_ids.size() );

	std::map< std::string, core::Real > retval;

	{
		core::Real ligand_rms_no_super= core::scoring::rmsd_no_super(
			test_pose,
			ref_pose,
			test_residue_ids,
			ref_residue_ids,
			core::scoring::is_ligand_heavyatom
		);

		std::string no_super_label("ligand_rms_no_super_");
		no_super_label += chain;
		if ( prefix != "" ) {
			no_super_label = prefix + "_" + no_super_label;
		}

		retval[ no_super_label ] = ligand_rms_no_super;
	}
	{
		core::Real ligand_rms_with_super= core::scoring::rmsd_with_super(
			test_pose,
			ref_pose,
			test_residue_ids,
			ref_residue_ids,
			core::scoring::is_ligand_heavyatom
		);

		std::string with_super_label("ligand_rms_with_super_");
		with_super_label += chain;
		if ( prefix != "" ) {
			with_super_label = prefix + "_" + with_super_label;
		}

		retval[ with_super_label ] = ligand_rms_with_super;
	}

	return retval;
}

} // namespace ligand_docking
} // namespace protocols
