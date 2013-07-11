// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/grafting/util.cc
/// @brief Utility functions for grafting
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)

//Unit headers
#include <protocols/grafting/util.hh>


//Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
//#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

//Basic headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.grafting");

namespace protocols {
namespace grafting {

using core::pose::Pose;

void
delete_region(Pose & pose, Size const start, Size const end){
	//Going to remove residues from pose. Resnum changes each time - so...
	Size const num_residues(end-start+1);
	TR << "Deleting "<< num_residues << " residues from " << start << " to "<< end << std::endl;

	for (core::Size i = 1; i <= num_residues; ++i) {
		//pose_.conformation().delete_polymer_residue(resnum);
		pose.conformation().delete_residue_slow(start);

	}
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
	for (core::Size resnum = start; resnum <= end; ++resnum){
		positions.push_back(resnum);
	}

	//Create the New foldtree to be applied in create_subpose
	//Simple FT for NOW.
	new_foldtree.simple_tree(length);


	core::pose::create_subpose(pose, positions, new_foldtree, piece);
	return piece;
}

Pose
replace_region(Pose const & to_pose, Pose const & from_pose, Size const from_pose_start_residue, Size const to_pose_start_residue, Size const insertion_length){
	Pose combined(to_pose);
	for (core::Size i= 0; i<= insertion_length-1; ++i) {
		core::Size pose_num = to_pose_start_residue+i;
		core::Size piece_num = from_pose_start_residue+i;
		core::conformation::Residue const & piece_rsd = from_pose.residue(piece_num);
		bool replace_backbone = true;

		combined.conformation().replace_residue(pose_num, piece_rsd, replace_backbone);

	}
	return combined;
}

///@author Steven Lewis smlewi@gmail.com
///@details brief inserts one pose into another pose, returning the product as a new value. The insert pose will be added immediately after insert_point and before insert_point_end.
Pose
insert_pose_into_pose(
											Pose const & scaffold_pose,
											Pose const & insert_pose,
											core::Size const insert_point,
											core::Size const insert_point_end
){
	//local copies
	core::pose::Pose scaffold(scaffold_pose);
	core::pose::Pose insert(insert_pose);
	
	//Get Disulfide information:
	utility::vector1< std::pair< core::Size, core::Size > > disulfide_pair_list;
	core::conformation::disulfide_bonds(insert.conformation(), disulfide_pair_list);
	
	core::kinematics::FoldTree original_scaffold_tree = scaffold.fold_tree();
	core::Size const insert_length = insert.total_residue();
	//strip termini variants from insert if necessary
	using core::pose::remove_variant_type_from_pose_residue;
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::LOWER_TERMINUS, 1);
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::UPPER_TERMINUS, insert_length);

	core::Size const insert_start(insert_point+1); //this will be first residue of the inflexible insert
	core::Size const insert_end(insert_point+insert_length); //this will be the last residue of the inflexible insert

	//Just in case someone wants to do something crazy and add something to termini on two different chains.
	//Not that it's been tested...
	core::pose::remove_variant_type_from_pose_residue(scaffold, core::chemical::LOWER_TERMINUS, insert_point);
	core::pose::remove_variant_type_from_pose_residue(scaffold, core::chemical::UPPER_TERMINUS, insert_point_end);
	
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
	for (core::Size i = 1; i <= insert_length; ++i){
		combined.insert_residue_by_jump(insert.residue(i), insert_point+i, insert_point);
	}

	//Set a fold tree.
	core::kinematics::FoldTree new_tree = core::kinematics::remodel_fold_tree_to_account_for_insertion(original_scaffold_tree, insert_point, insert_length);
	combined.fold_tree(new_tree);
	
	//Fix Disulfides.
	for (core::Size i = 1; i <= disulfide_pair_list.size(); ++i){
		core::Size new_resnum1 = disulfide_pair_list[i].first + insert_point;
		core::Size new_resnum2 = disulfide_pair_list[i].second + insert_point;

		//This is all we need, as the disulfides have already been detected in conformation. No need to mutate or optimize.
		combined.conformation().declare_chemical_bond(new_resnum1,"SG",new_resnum2,"SG");
	}
	
	return combined;
}



}
}
