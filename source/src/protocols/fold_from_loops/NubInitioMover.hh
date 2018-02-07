// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   NubInitioMover.hh
/// @brief  Ab Initio with a pre-folded segment (nub - pun intended)
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_NubInitioMover_hh
#define INCLUDED_protocols_fold_from_loops_NubInitioMover_hh

// Unit headers
#include <protocols/moves/Mover.hh>
#include <protocols/fold_from_loops/NubInitioMover.fwd.hh>
#include <protocols/fold_from_loops/utils/Nub.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueRanges.fwd.hh>
#include <core/io/silent/SilentFileData.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace fold_from_loops {

class NubInitioMover : public  protocols::moves::Mover {
public:

	/// @brief Empty Constructor
	NubInitioMover();
	/// @brief Destructor
	~NubInitioMover();

	/// @brief Main application of the NubInitioMover
	void apply( core::pose::Pose & pose ) override;

private:
	////////////////////////////////////////////////
	////////////////// MAIN STEPS //////////////////
	////////////////////////////////////////////////

	/// @brief Prepare the unfolded pose previous to folding
	void make_unfolded_pose( core::pose::Pose & pose );
	/// @brief Refold the unfolded pose
	core::Size refold_pose( core::pose::Pose & pose );
	/// @brief Finishup the folded structure
	void post_process( core::pose::Pose & pose );

private:
	////////////////////////////////////////////////
	/////////////// SUPPORT FUNCTIONS //////////////
	////////////////////////////////////////////////

	/// @brief Checks that input make sense for the particular run
	void sanity_check( core::pose::Pose const & pose );
	/// @brief List of Poses from the template that will become the UNFOLDED REGIONS.
	utility::vector1< core::pose::PoseOP > get_template_pieces( core::pose::Pose const & pose ) const;
	/// @brief Splits the template selected regions if needed.
	core::select::residue_selector::ResidueRangesOP
	make_template_ranges( core::select::residue_selector::ResidueRanges original ) const;
	// @brief Manages constraint requirements for the Pose.
	void manage_constraints( core::pose::Pose & pose );
	// @brief Run the ab initio simulation.
	core::scoring::ScoreFunctionOP apply_abinitio( core::pose::Pose & pose );
	/// @brief evaluate rmsd with template
	core::Real template_rmsd( core::pose::Pose & pose );
	/// @brief Counts how many residues in the pose DESIGN CHAIN are
	/// in contact with the inserted motif and labels them as CONTACT.
	core::Size count_contacts( core::pose::Pose & pose ) const;
	/// @brief Dumps the pose at centroid level in a silent file with the same file name
	/// and pose TAG but ending with "_CENTROID". Uses the provided score to evaluate it.
	void dump_centroid( core::pose::Pose const & pose, core::scoring::ScoreFunctionOP scorefxn );
	/// @brief fix the pose after refitting the sidechains
	void repack_minimize_disulfides( core::pose::Pose & pose );
	/// @brief Repack sidechains without alter bb
	void repack( core::pose::Pose & pose );

public:
	// -- SETTERS/GETTERS -- //
	/// @brief In rosettascripts, saves the user provided name.
	std::string prefix() const;
	void prefix( std::string const & prefix );
	/// @brief Store the unmodified pose as provided by apply.
	core::pose::Pose template_pose() const;
	core::select::residue_selector::ResidueSelectorCOP template_selector() const;
	void template_selector( core::select::residue_selector::ResidueSelectorCOP const & selector );
	/// @brief score function to use on repacking after going back to full atom
	core::scoring::ScoreFunctionOP fullatom_scorefxn() const;
	void fullatom_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );
	/// @brief object that contains and processes the unfolded pose info
	utils::NubOP nub() const;
	void nub( utils::NubOP nub );
	/// @brief when true, use constraint guided ab initio
	bool use_cst() const;
	void use_cst( bool pick );
	/// @brief when true, remove constraints transfered from the template into the inserted motif
	bool clear_motif_cst() const;
	void clear_motif_cst( bool pick );
	/// @brief limit of divergence allowed between the template structure and the ab initio output
	core::Real rmsd_threshold() const;
	void rmsd_threshold( core::Real threshold );
	/// @brief when true, include the residues of the insertion in the RMSD evaluation (do only when
	/// insertion sizes are the same)
	bool rmsd_include_motif() const;
	void rmsd_include_motif( bool pick );
	/// @brief when false, do not use residues without constraints in the RMSD evaluation
	bool rmsd_include_unconstrained() const;
	void rmsd_include_unconstrained( bool pick );
	/// @brief when true, try to repack disulfides
	bool repack_disulfides() const;
	void repack_disulfides( bool pick );
	/// @brief when true, allow bb movements to try to repack the disulfides
	bool disulfides_bb() const;
	void disulfides_bb( bool pick );
	/// @brief if disulfides_bb is true, how many residues around CYS are allowed to move?
	core::Size disulfides_side() const;
	void disulfides_side( core::Real value );
	/// @brief for rosettascripts, id of the fragmentsets in the datamap
	std::string fragments_id() const;
	void fragments_id( std::string const & name );
	/// @brief how many times do we try to fold?
	void max_trials( core::Size choice );
	core::Size max_trials() const;
	/// @brief weight for binder related score terms, to provide to the ab initio process
	void binder_weight( core::Real value );
	core::Real binder_weight() const;
	/// @brief weight for angle constraints, to provide to the ab initio process
	void angle_weight( core::Real value );
	core::Real angle_weight() const;
	/// @brief weight for dihedral constraints, to provide to the ab initio process
	void dihedral_weight( core::Real value );
	core::Real dihedral_weight() const;
	/// @brief value automatically applied to weight to improve alpha/beta fold
	void correction_weights( core::Real value );
	core::Real correction_weights() const;
	/// @brief if true, apply the correction weights
	bool use_correction_weights() const;
	void use_correction_weights( bool pick );
	/// @brief if true, generate silentfile with the centroid-level protein
	bool dump_centroid() const;
	void dump_centroid( bool pick );
	/// @brief if true, stop after generating the unfolded pose (for show and debug)
	bool drop_unfolded_pose() const;
	void drop_unfolded_pose( bool pick );
	/// @brief if true, design the designable residues (TEMPLATE & COLDSPOT labels)
	bool design() const;
	void design( bool pick );
	/// @brief if given, design all designable residues to the provided one (TEMPLATE & COLDSPOT labels)
	std::string residue_type() const;
	void residue_type( std::string pick );


public:
	// -- COMMON WORKING SELECTORS AND RESIDUERANGES -- //
	/// @brief Selects the chain considered as DESIGN by the NubInitioMover instance.
	core::select::residue_selector::ResidueSelectorOP design_chain_selector() const;
	/// @brief Selects ALL the chains NOT considered as DESIGN by the NubInitioMover instance.
	core::select::residue_selector::ResidueSelectorOP context_chain_selector() const;
	/// @brief INVERTS the INSERTION selector in the REFERENCE TEMPLATE POSE.
	core::select::residue_selector::ResidueSelectorOP used_template_selector() const;
	/// @brief Selects the residues that are MOTIF from the DESIGN CHAIN.
	core::select::residue_selector::ResidueSelectorOP work_motif_selector() const;
	/// @brief Selects the residues that are TEMPLATE from the DESIGN CHAIN.
	core::select::residue_selector::ResidueSelectorOP work_template_selector() const;
	/// @brief Selects the residues that are TEMPLATE or FLEXIBLE from the DESIGN CHAIN.
	core::select::residue_selector::ResidueSelectorOP bb_movable_selector() const;
	/// @brief Selects the residues that are TEMPLATE or COLDSPOT from the DESIGN CHAIN.
	core::select::residue_selector::ResidueSelectorOP chi_movable_selector() const;
	/// @brief Selects residues that are DISULFIDIZE from the DESIGN CHAIN.
	core::select::residue_selector::ResidueSelectorOP cys_design_selector() const;
	/// @brief Range of residues representing the INSERTION regions (checks that there is one)
	core::select::residue_selector::ResidueRanges template_insertion_ranges() const;
	/// @brief Range of residues kept from the TEMPLATE (checks that there is one)
	/// If terminals, adds N and C-terminal empty ranges (0, 0) if there is a N or C-terminal insertion
	core::select::residue_selector::ResidueRanges used_template_ranges( bool terminals = false ) const;

private:
	// -- DEFAULTS -- //
	static std::string default_prefix();
	static core::select::residue_selector::ResidueSelectorCOP default_template_selector();
	static core::scoring::ScoreFunctionOP default_fullatom_scorefxn();
	static bool default_use_cst();
	static bool default_clear_motif_cst();
	static core::Real default_rmsd_threshold();
	static bool default_rmsd_include_motif();
	static bool default_rmsd_include_unconstrained();
	static bool default_repack_disulfides();
	static bool default_disulfides_bb();
	static core::Size default_disulfides_side();
	static core::Size default_max_trials();
	static core::Real default_mc_binder_weight();
	static core::Real default_mc_angle_weight();
	static core::Real default_mc_dihedral_weight();
	static core::Real default_mc_correction_weights();
	static bool default_use_correction_weights();
	static bool default_dump_centroid();
	static bool default_drop_unfolded_pose();
	static bool default_design();
	static std::string default_residue_type();

public:
	// -- ROSETTASCRIPTS -- //
	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & reference_pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP clone() const override;
	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	// -- ATTRIBUTES -- //
	std::string prefix_;
	core::pose::Pose template_pose_;
	core::select::residue_selector::ResidueSelectorCOP template_selector_;
	core::scoring::ScoreFunctionOP fullatom_scorefxn_;
	core::scoring::ScoreFunctionOP abinitio_score_; // internal
	utils::NubOP nub_;
	bool use_cst_;
	bool clear_motif_cst_;
	core::Real rmsd_threshold_;
	bool rmsd_include_motif_;
	bool rmsd_include_unconstrained_;
	bool repack_disulfides_;
	bool disulfides_bb_;
	core::Size disulfides_side_;
	std::string fragments_id_;
	core::Real trials_;
	core::Real max_trials_;
	core::Real mc_binder_weight_;
	core::Real mc_angle_weight_;
	core::Real mc_dihedral_weight_;
	core::Real mc_correction_weights_;
	bool use_correction_weights_;
	bool has_betas_;  // internal
	bool has_alphas_; // internal
	bool dump_centroid_;
	bool drop_unfolded_pose_;
	bool design_;
	std::string residue_type_;
	core::io::silent::SilentFileDataOP silent_score_file_;

};

}
}

#endif
