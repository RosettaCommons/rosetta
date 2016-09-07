// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/grafting/util.cc
/// @brief Utility functions for grafting
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)

//Unit headers
#include <protocols/grafting/util.hh>


//Project headers
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

//Protocols
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/util.hh>
#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>

//#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

//Basic headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.grafting.util" );

namespace protocols {
namespace grafting {

using namespace core;

using namespace core::id;
using core::pose::Pose;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
using core::kinematics::MoveMapCOP;
using core::Size;
using protocols::loops::Loops;
using core::scoring::ScoreFunctionCOP;

void
delete_region(Pose & pose, Size const start, Size const end){
	//Going to remove residues from pose. Resnum changes each time - so...
	Size const num_residues(end-start+1);
	TR << "Deleting "<< num_residues << " residues from " << start << " to "<< end << std::endl;

	if ( num_residues == 0 ) {
		TR <<"Nothing to delete..." << std::endl;
		return;
	}
	pose.delete_residue_range_slow(start, end);
	pose.conformation().detect_disulfides();
	//for (core::Size i = 1; i <= num_residues; ++i) {
	//pose_.conformation().delete_polymer_residue(resnum);
	//pose.conformation().delete_residue_slow(start);

	//}
	if ( pose.pdb_info() ) pose.pdb_info()->obsolete(false);
}


Pose
return_region(Pose & pose, Size const start, Size const end){

	//Copy residues into a new pose.  Uses create_subpose.
	Pose piece;
	core::kinematics::FoldTree new_foldtree;
	new_foldtree.clear();
	Size const length = end - start+1;
	TR << "Returning "<< length <<" residues from "<< start <<" to " << end <<std::endl;

	//Get positions.
	utility::vector1<Size> positions;
	for ( core::Size resnum = start; resnum <= end; ++resnum ) {
		positions.push_back(resnum);
	}

	//Create the New foldtree to be applied in create_subpose
	//Simple FT for NOW.
	new_foldtree.simple_tree(length);


	core::pose::create_subpose(pose, positions, new_foldtree, piece);

	//Create subpose results in a NULL PDB_Info.  We now need a new one.
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo(piece.total_residue()) );
	piece.pdb_info(pdb_info);

	piece.pdb_info()->copy(*(pose.pdb_info()), start, end, 1); //Should be an option directly within subpose
	piece.pdb_info()->obsolete(false);
	piece.conformation().detect_disulfides();

	return piece;
}

Pose
replace_region(
	Pose const & from_pose,
	Pose const & to_pose,
	Size const from_pose_start_residue,
	Size const to_pose_start_residue,
	Size const insertion_length,
	bool copy_pdbinfo /*false*/)
{
	Pose combined(to_pose);
	for ( core::Size i= 0; i<= insertion_length-1; ++i ) {
		core::Size pose_num = to_pose_start_residue+i;
		core::Size piece_num = from_pose_start_residue+i;
		core::conformation::Residue const & piece_rsd = from_pose.residue(piece_num);
		bool replace_backbone = true;

		combined.conformation().replace_residue(pose_num, piece_rsd, replace_backbone);

	}


	if ( copy_pdbinfo && from_pose.pdb_info() && to_pose.pdb_info() ) {
		combined.pdb_info()->copy(*(from_pose.pdb_info()), from_pose_start_residue, from_pose_start_residue + insertion_length -1, to_pose_start_residue);
		combined.pdb_info()->obsolete(false);
	}
	combined.conformation().detect_disulfides();
	return combined;
}

Pose
insert_pose_into_pose(
	Pose const & scaffold_pose,
	Pose const & insert_pose,
	core::Size const insert_point,
	bool copy_pdbinfo /*false*/

){
	return insert_pose_into_pose(scaffold_pose, insert_pose, insert_point, insert_point+1, copy_pdbinfo);
}

/// @author Steven Lewis smlewi@gmail.com
/// @details brief inserts one pose into another pose, returning the product as a new value. The insert pose will be added immediately after insert_point and before insert_point_end.
Pose
insert_pose_into_pose(
	Pose const & scaffold_pose,
	Pose const & insert_pose,
	core::Size const insert_point,
	core::Size const insert_point_end,
	bool copy_pdbinfo /*false*/
){
	//local copies
	core::pose::Pose scaffold(scaffold_pose);
	core::pose::Pose insert(insert_pose);

	//TR << "Pre Insertion: " << std::endl;
	//scaffold_pose.constraint_set()->show(TR);
	//scaffold.constraint_set()->show(TR);

	//Get Disulfide information:
	utility::vector1< std::pair< core::Size, core::Size > > disulfide_pair_list;
	core::conformation::disulfide_bonds(insert.conformation(), disulfide_pair_list);

	core::kinematics::FoldTree original_scaffold_tree = scaffold.fold_tree();
	core::Size const insert_length = insert.total_residue();
	//strip termini variants from insert if necessary
	using core::pose::remove_variant_type_from_pose_residue;
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::LOWER_TERMINUS_VARIANT, 1);
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::UPPER_TERMINUS_VARIANT, insert_length);

	core::Size const insert_start(insert_point+1); //this will be first residue of the inflexible insert
	core::Size const insert_end(insert_point+insert_length); //this will be the last residue of the inflexible insert

	//Just in case someone wants to do something crazy and add something to termini on two different chains.
	//Not that it's been tested...
	core::pose::remove_variant_type_from_pose_residue(scaffold, core::chemical::LOWER_TERMINUS_VARIANT, insert_point);
	core::pose::remove_variant_type_from_pose_residue(scaffold, core::chemical::UPPER_TERMINUS_VARIANT, insert_point_end);

	TR << "insert_point " << insert_point << std::endl;
	TR << "insert_start " << insert_start << std::endl;
	TR << "insert_end " << insert_end << std::endl;
	TR << "insert_length " << insert_length << std::endl;

	//Fold tree allows insertion into scaffold (all via jump)
	using core::kinematics::Edge;
	core::pose::Pose combined(scaffold);
	core::kinematics::FoldTree inserting_tree(combined.total_residue());
	inserting_tree.clear();
	inserting_tree.add_edge(Edge(1, insert_point, Edge::PEPTIDE));
	inserting_tree.add_edge(Edge(insert_point_end, combined.total_residue(), Edge::PEPTIDE));
	inserting_tree.add_edge(Edge(insert_point, insert_point_end, 1));
	inserting_tree.reorder(1);
	TR << inserting_tree << std::endl;
	combined.fold_tree(inserting_tree);

	//insert insert residues by jump (just to get them in)
	for ( core::Size i = 1; i <= insert_length; ++i ) {
		combined.insert_residue_by_jump(insert.residue(i), insert_point+i, insert_point);
		//TR << "Added residue: "<< i << std::endl;
		//combined.constraint_set()->show(TR);
		//TR << std::endl;
	}

	//Set a fold tree.
	core::kinematics::FoldTree new_tree = core::kinematics::remodel_fold_tree_to_account_for_insertion(original_scaffold_tree, insert_point, insert_length);
	combined.fold_tree(new_tree);

	if ( copy_pdbinfo && insert_pose.pdb_info() && scaffold_pose.pdb_info() ) {
		//TR << "Start PDBInfo size: " << combined.pdb_info()->nres() << std::endl;
		//TR << "Insert PDBInfo Size: " << insert_pose.pdb_info()->nres() << std::endl;
		combined.pdb_info()->copy(*(insert_pose.pdb_info()), 1, insert_pose.total_residue(), insert_point+1);
		combined.pdb_info()->obsolete(false);
		//TR << "Final PDBInfo Size: " << combined.pdb_info()->nres() << std::endl;
	}

	//Fix Disulfides.
	for ( core::Size i = 1; i <= disulfide_pair_list.size(); ++i ) {
		core::Size new_resnum1 = disulfide_pair_list[i].first + insert_point;
		core::Size new_resnum2 = disulfide_pair_list[i].second + insert_point;

		//This is all we need, as the disulfides have already been detected in conformation. No need to mutate or optimize.
		combined.conformation().declare_chemical_bond(new_resnum1,"SG",new_resnum2,"SG");
	}

	//TR << "Post insertion: " << std::endl;
	//combined.constraint_set()->show(TR);

	return combined;
}

//////////////////////////////////////////////////////
void
repack_connection_and_residues_in_movemap(
	core::pose::Pose & pose, core::scoring::ScoreFunctionCOP fa_scorefxn,
	core::Size const start, core::Size const end, core::kinematics::MoveMapCOP movemap){

	using protocols::toolbox::task_operations::RestrictToMoveMapChiOperationOP;
	using namespace core::pack::task::operation;
	using namespace core::pack::task;
	MoveMapOP local =movemap->clone();
	local->set_chi(start, true);
	local->set_chi(start+1, true);
	local->set_chi(end, true);
	local->set_chi(end-1, true);

	RestrictToMoveMapChiOperationOP mm_task_op( new protocols::toolbox::task_operations::RestrictToMoveMapChiOperation(local) );
	mm_task_op->set_include_neighbors(false);

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ));
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ));
	tf->push_back(mm_task_op);

	PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover(fa_scorefxn, task) );
	packer->apply(pose);
}

void
repack_connection_and_residues_in_movemap_and_piece(
	core::pose::Pose & pose, core::scoring::ScoreFunctionCOP fa_scorefxn,
	core::Size const start, core::Size const end, core::kinematics::MoveMapCOP movemap){

	using protocols::toolbox::task_operations::RestrictToMoveMapChiOperationOP;
	using namespace core::pack::task::operation;
	using namespace core::pack::task;

	MoveMapOP local = movemap->clone();
	local->set_chi(start, true);
	local->set_chi(start+1, true);
	local->set_chi(end, true);
	local->set_chi(end-1, true);
	for ( Size i = start+2; i <= end-2; ++i ) {
		local->set_chi(i, true);
	}

	RestrictToMoveMapChiOperationOP mm_task_op( new protocols::toolbox::task_operations::RestrictToMoveMapChiOperation(local) );
	mm_task_op->set_include_neighbors(false);

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ));
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ));
	tf->push_back(mm_task_op);

	PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover(fa_scorefxn, task) );
	packer->apply(pose);
}

void
repack_connection_and_residues_in_movemap_and_piece_and_neighbors(
	core::pose::Pose & pose, ScoreFunctionCOP fa_scorefxn,
	core::Size const start, core::Size const end, MoveMapCOP movemap, core::Real neighbor_dis)
{
	using protocols::toolbox::task_operations::RestrictToMoveMapChiOperationOP;
	using namespace core::pack::task::operation;
	using namespace core::pack::task;

	MoveMapOP local = movemap->clone();
	local->set_chi(start, true);
	local->set_chi(start+1, true);
	local->set_chi(end, true);
	local->set_chi(end-1, true);
	for ( Size i = start+2; i <= end-2; ++i ) {
		local->set_chi(i, true);
	}

	RestrictToMoveMapChiOperationOP mm_task_op( new protocols::toolbox::task_operations::RestrictToMoveMapChiOperation(local) );
	mm_task_op->set_include_neighbors(true);
	mm_task_op->set_cutoff_distance(neighbor_dis);

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ));
	tf->push_back(TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ));
	tf->push_back(mm_task_op);

	PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);
	protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover(fa_scorefxn, task) );
	packer->apply(pose);
}

void
superimpose_overhangs_heavy(Pose const & pose, Pose & piece, bool ca_only, Size start, Size end, Size Nter_overhang, Size Cter_overhang){

	if ( Nter_overhang == 0 && Cter_overhang ==0 ) {
		return;
	}
	TR <<"Superimposing overhang residues" << std::endl;
	AtomID_Map < AtomID > atoms_to_superimpose;
	initialize_atomid_map( atoms_to_superimpose, piece, core::id::BOGUS_ATOM_ID );
	//Remove termini
	remove_lower_terminus_type_from_pose_residue(piece, 1);
	remove_upper_terminus_type_from_pose_residue(piece, piece.total_residue());

	//Nter residues/atoms
	for ( core::Size piece_res_num=1; piece_res_num<=Nter_overhang; ++piece_res_num ) {
		core::Size pose_res_num = start - Nter_overhang + piece_res_num;
		if ( ca_only ) {
			AtomID const atom_piece(piece.residue(piece_res_num).atom_index("CA"), piece_res_num);
			AtomID const atom_pose(pose.residue(pose_res_num).atom_index("CA"), pose_res_num);
			atoms_to_superimpose[ atom_piece ]=atom_pose;
		} else {
			for ( core::Size atom_num=1; atom_num <=4; atom_num++ ) {
				AtomID const atom_piece(atom_num, piece_res_num);
				AtomID const atom_pose(atom_num, pose_res_num);
				atoms_to_superimpose[ atom_piece ]=atom_pose;
			}
		}
	}

	//Cter residues/atoms
	core::Size i = 0;
	for ( core::Size piece_res_num=(piece.total_residue() - Cter_overhang + 1); piece_res_num<=piece.total_residue(); ++piece_res_num ) {
		core::Size pose_res_num = end + i;
		if ( ca_only ) {
			AtomID const atom_piece(piece.residue(piece_res_num).atom_index("CA"), piece_res_num);
			AtomID const atom_pose(pose.residue(pose_res_num).atom_index("CA"), pose_res_num);
			atoms_to_superimpose[ atom_piece ]=atom_pose;
		} else {
			for ( core::Size atom_num=1; atom_num <=4; ++atom_num ) {
				AtomID const atom_piece(atom_num, piece_res_num);
				AtomID const atom_pose( atom_num, pose_res_num);
				atoms_to_superimpose[ atom_piece ]=atom_pose;
			}
		}
		++i;
	}

	core::Real rms(core::scoring::superimpose_pose(piece, pose, atoms_to_superimpose));
	TR.Debug << "RMS: " << rms << std::endl;
	return;
}

void
delete_overhang_residues(Pose & piece, Size Nter_overhang, Size Cter_overhang){

	//Nter
	for ( core::Size i=1; i<=Nter_overhang; ++i ) {
		//piece_.conformation().delete_polymer_residue(piece_res_num);
		piece.conformation().delete_residue_slow(1);
	}

	//Cter
	core::Size res_num_start = (piece.total_residue() - Cter_overhang + 1);
	for ( core::Size i=1; i<=Cter_overhang; ++i ) {
		//piece_.conformation().delete_polymer_residue(piece_res_num);
		piece.conformation().delete_residue_slow(res_num_start);
	}
	piece.pdb_info()->obsolete(false);

}


MoveMapOP
combine_movemaps_post_insertion(MoveMapCOP scaffold_mm, MoveMapCOP insert_mm,
	Size start, Size original_end,
	Size insertion_length, Size cter_overhang_start)
{
	using namespace core::kinematics;
	MoveMapOP mm( new MoveMap() );
	TR<<"Combining movemaps"<<std::endl;
	//typedef core::id::TorsionID TorsionID;
	typedef std::pair< Size, core::id::TorsionType > MoveMapTorsionID;

	//Assert that a piece is set.

	//Copy Nterminal MoveMapTorsionIDs
	for ( auto it=scaffold_mm->movemap_torsion_id_begin(), it_end=scaffold_mm->movemap_torsion_id_end(); it !=it_end; ++it ) {
		//Scaffold to new MM
		if ( it->first.first<=start ) {
			MoveMapTorsionID new_id = MoveMapTorsionID(it->first.first, it->first.second);
			mm->set(new_id, it->second);
		}
	}

	//Set insert residues
	for ( auto it=insert_mm->movemap_torsion_id_begin(), it_end=insert_mm->movemap_torsion_id_end(); it !=it_end; ++it ) {
		//Add from start_
		Size new_resnum = start + it->first.first -  cter_overhang_start;
		MoveMapTorsionID new_id = MoveMapTorsionID(new_resnum, it->first.second);
		mm->set(new_id, it->second);
	}

	//Set Cterminal residues after insert. We may have a deletion then an insertion that we need to change the numbers for.
	Size deleted_residues = original_end-start-1;
	for ( auto it=scaffold_mm->movemap_torsion_id_begin(), it_end=scaffold_mm->movemap_torsion_id_end(); it !=it_end; ++it ) {
		//Check if residue exists in movemap.  Copy that info no matter what it is to the new movemap.
		if ( it->first.first >=original_end ) {
			//Add insertion length
			Size new_resnum = it->first.first - deleted_residues+insertion_length;
			MoveMapTorsionID new_id = MoveMapTorsionID(new_resnum, it->first.second);
			mm->set(new_id, it->second);
		}
	}
	mm->show(TR.Debug);
	return mm;
}

core::Real
perturb_backbone_for_test(Pose& pose, MoveMapOP mm){
	Pose pre_test(pose);
	TR <<"Randomizing residues in movemap." <<std::endl;
	protocols::simple_moves::SmallMover randomize(mm, 10000, 500); //Huge KT to randomize the residues
	randomize.angle_max( 'H', 180.0 );
	randomize.angle_max( 'E', 180.0 );
	randomize.angle_max( 'L', 180.0 );
	randomize.apply(pose);

	Real RMS = core::scoring::bb_rmsd_including_O(pre_test, pose);
	TR <<"All BB Heavy Atom RMSD "<< RMS <<std::endl;
	return RMS;
}

void
add_cutpoint_variants_for_ccd(Pose & pose, Loops const & loops){
	core::Size i = 0;
	 for ( auto const & loop : loops ) {
		++i;
		TR <<"Loop: "<< i << std::endl;
		TR << "Add variant to: " << loop.cut() << std::endl;
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_LOWER, loop.cut());
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_UPPER, loop.cut()+1);
	}

}

void
remove_cutpoint_variants_for_ccd(Pose & pose, Loops const & loops){
	 for ( auto const & loop : loops ) {
		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_LOWER, loop.cut());
		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_UPPER, loop.cut()+1);
	}
}

bool
graft_closed(Pose & pose, Loops & loops){
	 for ( auto const & loop : loops ) {
		std::pair<bool, Size> peptide_bond_issues = protocols::loops::has_severe_pep_bond_geom_issues(pose, loop.cut(), true, true, 1.5, 15, 15);
		if ( peptide_bond_issues.first == true ) {
			return false;
		}
	}

	return true;
}

void
idealize_combined_pose(Pose & combined, MoveMapOP movemap, Size start, Size insert_start, Size insert_end, Size Nter_loop_start, Size Cter_loop_end, bool idealize_insert /*false*/){
	///////////////////////////////////////Idealize////////////////////////////////////////////////////////
	//this code also resets conformation variables: omegas to 180, newly made connections phi or psi to reasonable
	//edges of insert will be somewhat mobile inside minimization (small and CCD moves will ignore it)

	//Order of idealization matters here.
	//Nter regions
	for ( core::Size i = Nter_loop_start; i<=start; ++i ) {
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		if ( movemap->get_bb(i) ) {
			combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
			combined.set_omega(i, 180);
			TR << "idealized  " << i << std::endl;
		}
	}
	TR << "ideal " << insert_start << std::endl;

	//Set individual torsions ON in the movemap for the start and end of the insert

	combined.set_phi(insert_start, -60);

	//insert regions
	if ( idealize_insert ) {
		for ( core::Size i = start +1; i<=insert_end; ++i ) {
			//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
			if ( movemap->get_bb(i) ) {
				combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
				combined.set_omega(i, 180);
				TR << "idealized " << i << std::endl;
			}
		}
	}


	combined.conformation().insert_ideal_geometry_at_polymer_bond(insert_end);

	TR << "ideal " << insert_end << std::endl;
	combined.set_omega(insert_end, 180);//Note, upper Nter idealize loop takes care of omega for insert_start-1
	combined.set_psi(insert_end, -40);

	//Cter regions

	for ( core::Size i=insert_end+1; i<=Cter_loop_end; ++i ) {
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		if ( movemap->get_bb(i) ) {
			combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
			combined.set_omega(i, 180);
			TR << "idealized " << i << std::endl;
		}
	}

}
/////////////////////////////////////////////////////////////////////////////////////////////////
///FOLDTREE SETUP.   options depending on how you want your graft algorithm to work!
///---Indicates Flexible regions, | indicates cutpoint. Arrows are direction of ARMs used to close the loop in conjunction with algorithm. (CCD, KIC)
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece----> | Cter_loop_end****
/// Insert will move in cartesian space
/// @params lower_cutpoint for CCD and loops is Cter_loop_end-1
///
void
setup_single_loop_single_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, bool loop_modeling){
	using core::kinematics::Edge;

	//Setup the offset so the anchor of the loop is correct for either CCD loop closure or vanilla loop modeling.
	Size anchor_offset = 1;
	if ( loop_modeling ) { anchor_offset=2;}

	core::Size const cutpoint_lower(Cter_loop_end-1); //cutpoint at the end of the loop
	core::Size const cutpoint_upper(Cter_loop_end); //cutpoint at the end of the loop
	core::Size const loop_start_foldtree_anchor(Nter_loop_start-anchor_offset); //this is the N-terminal jump anchor for the loop
	core::Size const loop_end_foldtree_anchor(Cter_loop_end+anchor_offset); //C-terminal jump anchor

	TR.Debug << "loop_start_foldtree_anchor " << loop_start_foldtree_anchor << std::endl;
	TR.Debug << "cutpoint_lower " << cutpoint_lower << std::endl;
	TR.Debug << "cutpoint_upper " << cutpoint_upper << std::endl;
	TR.Debug << "loop_end_foldtree_anchor " << loop_end_foldtree_anchor << std::endl;

	core::kinematics::FoldTree remodeling_tree(pose.total_residue());
	remodeling_tree.clear();
	remodeling_tree.add_edge(Edge(1, loop_start_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, cutpoint_lower, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(cutpoint_upper, loop_end_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_end_foldtree_anchor, pose.total_residue(), Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, loop_end_foldtree_anchor, 1));
	remodeling_tree.reorder(1);
	TR << remodeling_tree << std::endl;
	pose.fold_tree(remodeling_tree);
}


//////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece | <----Nter_loop_end****
/// Insert will move in cartesian space
/// @params lower_cutpoint for CCD and loops is end_-1
///
void
setup_single_loop_double_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, Size end, bool loop_modeling){
	//Note: Differs from single arm by moving the cutpoint.
	TR <<"Setting up single loop, double arm remodeling foldtree"<<std::endl;
	using core::kinematics::Edge;

	//Setup the offset so the anchor of the loop is correct for either CCD loop closure or vanilla loop modeling.
	Size anchor_offset = 1;
	if ( loop_modeling ) { anchor_offset=2;}

	core::Size const cutpoint_lower(end-1);//Cutpoint is at the end of the insert
	core::Size const cutpoint_upper(end);
	core::Size const loop_start_foldtree_anchor(Nter_loop_start-anchor_offset); //this is the N-terminal jump anchor for the loop
	core::Size const loop_end_foldtree_anchor(Cter_loop_end+anchor_offset); //C-terminal jump anchor

	TR.Debug << "loop_start_foldtree_anchor " << loop_start_foldtree_anchor << std::endl;
	TR.Debug << "cutpoint_lower " << cutpoint_lower << std::endl;
	TR.Debug << "cutpoint_upper " << cutpoint_upper << std::endl;
	TR.Debug << "loop_end_foldtree_anchor " << loop_end_foldtree_anchor << std::endl;

	core::kinematics::FoldTree remodeling_tree(pose.total_residue());
	remodeling_tree.clear();
	remodeling_tree.add_edge(Edge(1, loop_start_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, cutpoint_lower, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(cutpoint_upper, loop_end_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_end_foldtree_anchor, pose.total_residue(), Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, loop_end_foldtree_anchor, 1));
	remodeling_tree.reorder(1);
	TR << remodeling_tree << std::endl;
	pose.fold_tree(remodeling_tree);
}


} //grafting
} //protocols
