// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh
/// @details Holds flags, options, and pose variables for the Glycosylation Protocol. Convenient object
/// to pass glycopeptide_docking options around.
/// @author Yashes Srinivasan (yashess@gmail.com)


#ifndef INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingFlags_hh
#define INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingFlags_hh

#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace glycopeptide_docking {

///@brief Holds flags, options, and pose variables for the Glycosylation Protocol.
class GlycopeptideDockingFlags : public utility::VirtualBase {

public:

	GlycopeptideDockingFlags();

	GlycopeptideDockingFlagsOP clone() const;

	/// @brief Get index (internal numbering) of the first residue of the substrate
	core::Size first_residue_substrate() const;

	/// @brief Get index (internal numbering) of the last residue of the substrate
	core::Size last_residue_substrate() const;

	/// @brief Get index (internal numbering) of the
	/// fold tree anchor (fold tree goes out in all direction from this residue).
	/// Only important when fold tree is specified as 'outward'.
	core::Size anchor_residue_substrate() const;

	/// @details Get index (internal numbering) of the site of glycosylation
	/// under investigation. (There may be other glycosylated residues on the peptide)
	core::Size glycosylation_residue_substrate() const;

	/// @brief Get jump number of the substrate
	core::Size jump_num_substrate() const;

	/// @brief Get first residue of the enzyme
	core::Size first_residue_enzyme() const;

	/// @brief Get last residue of the enzyme
	core::Size last_residue_enzyme() const;

	/// @brief Set variables from specified command line options.
	void
	setup_from_options();

	/// @brief Set enzyme chain.
	void
	set_enzyme_chain( core::pose::Pose const &pose );

	/// @brief Set anchor residue. Residue from which the foldtree
	/// emnates outwards when fold tree type is set as "outward"
	void
	set_anchor_residue(core::pose::Pose const &pose,std::string const &special_residue);

	/// @brief Set anchor residue with internal numbering
	void
	set_anchor_residue(core::Size in){
		anchor_residue_substrate_ = in;
	}

	/// @brief Get anchor residue with internal numbering
	core::Size
	get_anchor_residue(){
		return anchor_residue_substrate_;
	}

	/// @brief Set glycosylation residue from pose and pdb numbering
	void
	set_glycosylation_residue(core::pose::Pose const &pose,std::string const &special_residue);

	/// @brief Set glycosylation residue with internal numbering
	void
	set_glycosylation_residue(core::Size in){
		residue_to_glycosylate_ = in;
	}

	/// @brief Get glycosylation residue with internal numbering
	core::Size
	get_glycosylation_residue(){
		return residue_to_glycosylate_;
	}

	/*void
	set_residue_to_randomize(core::pose::Pose const &pose,std::string special_residue);
	*/

	/// @brief Set substrate chain.
	void
	set_substrate_chain( core::pose::Pose const &pose );

	/// @brief Show set options.
	void
	show() const; // For debugging purposes

	/// @brief Set whether to run in score_only mode. Protocol is not run.
	/// Only score is evaluated with distance, rmsd and interaction energy metrics.
	void set_score_only(bool in){
		score_only_ = in;
	}

	/// @brief Get whether to run in sampling or score-only mode.
	bool score_only(){
		return score_only_;
	}

	/// @brief Get whether to output interim debug pdbs during low and high resolution
	/// cycles.
	bool debug_pdbs(){
		return debug_pdbs_;
	}

	/// @brief Set whether to output interim debug pdbs during low and high resolution
	/// cycles.
	void set_debug_pdbs(bool in){
		debug_pdbs_=in;
	}

	/// @brief Set interface distance
	void set_interface_distance(core::Real distance){
		interface_distance_ = distance;
	}

	/// @brief Get interface distance
	core::Real get_interface_distance(){
		return interface_distance_;
	}

	/// @brief Set low and high res refinement options
	void set_glycosylation_refinement(bool low, bool high){
		glycosylation_high_res_refinement_=low;
		glycosylation_high_res_refinement_=high;
	}

	/// @brief Get whether low refinement is on
	bool low_res_refinement(){
		return glycosylation_low_res_refinement_;
	}

	/// @brief Get whether high refinement is on
	bool high_res_refinement(){
		return glycosylation_high_res_refinement_;
	}

	/// @brief Get how many outer MC cycles to run for high-res refinement
	core::Size high_res_outer_cycles(){
		return high_res_outer_cycles_;
	}

	/// @brief Get how many outer MC cycles to run for low-res refinement
	core::Size low_res_outer_cycles(){
		return low_res_outer_cycles_;
	}

	/// @brief Get how many inner MC cycles to run for low-res refinement
	core::Size low_res_inner_cycles(){
		return low_res_inner_cycles_;
	}

	/// @brief Get residue index of sugar donor
	core::Size get_sugar_donor(){
		return donor_;
	}

	/// @brief Set residue index of sugar donor
	void set_sugar_donor(core::Size in){
		donor_=in;
	}

	/// @brief Get substrate type: Only 'peptide' type is supported for now
	std::string substrate_type(){
		return substrate_type_;
	}

	/// @brief Set substrate type: Only 'peptide' type is supported for now
	/// TODO: lipids - if there are interesting applications
	void set_substrate_type(std::string const &in){
		substrate_type_=in;
	}

	//void set_idealize_rings(bool in){
	//idealize_rings_=in;
	//}

	/// @brief Get upstream chain for 'outward' fold tree
	std::string get_upstream_chain(){
		return upstream_chains_;
	}

	/// @brief Get downstream chain for 'outward' fold tree
	std::string get_downstream_chain(){
		return downstream_chains_;
	}

	/// @brief Set upstream chain for 'outward' fold tree
	void set_upstream_chain(std::string const &chain){
		upstream_chains_=chain;
	}

	/// @brief Set downstream chain for 'outward' fold tree
	void set_downstream_chain(std::string const &chain){
		downstream_chains_=chain;
	}

	/// @brief Set fold tree type: outward OR docking
	void set_fold_tree_type(std::string const &in){
		tree_type_=in;
	}

	/// @brief Get fold tree type
	std::string tree_type(){
		return tree_type_;
	}

	/// @brief Set how frequently the interface must be sampled:
	/// Eg. every 3rd cycle
	void set_interface_moves(core::Size n){
		nevery_interface_moves_=n;
	}

	/// @brief Get how frequently the interface is sampled
	core::Size get_interface_moves(){
		return nevery_interface_moves_;
	}

	/// @brief Set cycles for backbone moves in high resolution
	void set_backbone_moves(core::Size n){
		ntotal_backbone_moves_=n;
	}

	/// @brief Get cycles for backbone moves in high resolution
	core::Size get_backbone_moves(){
		return ntotal_backbone_moves_;
	}

	/// @brief On anchor residue sidechain repacking on or off
	void set_anchor_repacking(bool in){
		prevent_anchor_repacking_=in;
	}

	/// @brief Get whether anchor repacking on or off
	bool pack_anchor(){
		return prevent_anchor_repacking_;
	}

	/// @brief Set whether to sample glycan torsional dofs with glycan sampler
	void set_allow_glycan_torsion_moves(bool in){
		allow_glycan_torsion_moves_=in;
	}

	/// @brief Get whether to sample glycan torsional dofs with glycan sampler
	bool get_allow_glycan_torsion_moves(){
		return allow_glycan_torsion_moves_;
	}

	/// @brief Set name of constraints file
	/// Useful for setting constraints to preserve any catalytic motifs
	void set_constraints_file(std::string const & filename){
		constraints_=filename;
	}

	/// @brief Get name of constraints file
	std::string get_constraints_file(){
		return constraints_;
	}

	/// @brief Set whether to randomize initial torsions
	void set_randomize_substrate_torsions(bool in){
		randomize_substrate_torsions_=in;
	}

	/// @brief Get whether to randomize initial torsions
	bool randomize_substrate_torsions(){
		return randomize_substrate_torsions_;
	}

	/* void set_randomize_residue(core::Size in){
	residue_to_randomize_=in;
	}

	void enable_randomize_residue(bool in){
	randomize_residue_=in;
	}

	bool randomize_residue(){
	return randomize_residue_;
	}

	core::Size get_randomize_residue(){
	return residue_to_randomize_;
	}*/

	/// @brief Get whether to sample substrate backbone
	/// TODO: Add similar option for enzyme
	bool enable_backbone_moves_substrate(){
		return enable_backbone_moves_pp_;
	}

	/// @brief Set whether to sample substrate backbone
	void set_backbone_moves_substrate(bool in){
		enable_backbone_moves_pp_=in;
	}

	/// @brief Get whether to output additional distance metrics to scorefile
	/// - Get the same output as publication Mahajan et. al. 2020 Biorxiv
	bool additional_metrics(){
		return additional_metrics_;
	}

	/// @brief Enable writing of additional metrics to scorefile
	void enable_additional_metrics(bool in){
		additional_metrics_ = in;
	}

	/// @brief Set residue indices (int numbering) to pre-glycosylate before starting sampling
	void set_residues_to_preglycosylate(utility::vector1<core::Size> in){
		preglycosylate_residues_ = in;
	}

	/// @brief Get residue indices (int numbering) to pre-glycosylate before starting sampling
	utility::vector1<core::Size> get_residues_to_preglycosylate(){
		return preglycosylate_residues_;
	}

	/// @brief Set names of sugars for pre-glycosylate before starting sampling
	/// #sugar names must match #residue indices for preglycosylation.
	void set_sugars_to_preglycosylate(utility::vector1<std::string> in){
		preglycosylate_sugars_ = in;
	}

	/// @brief Get names of sugars for pre-glycosylate before starting sampling
	/// #sugar names must match #residue indices for preglycosylation.
	utility::vector1<std::string> get_sugars_to_preglycosylate(){
		return preglycosylate_sugars_;
	}



private:
	core::Size first_residue_substrate_;
	core::Size last_residue_substrate_;
	core::Size anchor_residue_substrate_;
	core::Size residue_to_glycosylate_;
	core::Size jump_num_substrate_;

	core::Size first_residue_enzyme_;
	core::Size last_residue_enzyme_;

	bool glycosylation_high_res_refinement_=true; //done
	bool glycosylation_low_res_refinement_=true;  //done
	std::string constraints_="";

	core::Size high_res_outer_cycles_=10;
	core::Size low_res_outer_cycles_=30;
	core::Size low_res_inner_cycles_=50;

	std::string substrate_type_="peptide"; //done
	std::string tree_type_="docking"; //done
	bool randomize_substrate_torsions_=false;
	bool enable_backbone_moves_pp_=true;
	std::string downstream_chains_=""; //done
	std::string upstream_chains_=""; //done
	bool prevent_anchor_repacking_=true;//done
	core::Size donor_=1;//done dummy value
	bool debug_pdbs_=false; //done
	core::Size nevery_interface_moves_=3;//done
	core::Size ntotal_backbone_moves_=30;//done
	bool allow_glycan_torsion_moves_=false;//done
	core::Real interface_distance_=8.0;//done
	bool score_only_=false; //done
	bool additional_metrics_=false;
	utility::vector1<core::Size> preglycosylate_residues_;
	utility::vector1<std::string> preglycosylate_sugars_;
	//core::Size residue_to_randomize_;
	//bool randomize_residue_=false;
};


} //glycopeptide_docking
} //protocols



#endif //INCLUDED_protocols_glycopeptide_docking_GlycopeptideDockingFlags_hh





