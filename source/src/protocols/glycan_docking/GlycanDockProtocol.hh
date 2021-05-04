// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/glycan_docking/GlycanDockProtocol.hh
/// @brief
/// @author Morgan Nance (morganlnance@gmail.com)

#ifndef INCLUDED_protocols_glycan_docking_GlycanDockProtocol_HH
#define INCLUDED_protocols_glycan_docking_GlycanDockProtocol_HH

// Unit headers
#include <protocols/glycan_docking/GlycanDockProtocol.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.fwd.hh>

// Protocol headers
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>
#include <protocols/carbohydrates/LinkageConformerMover.fwd.hh>
#include <protocols/constraint_movers/ConstraintSetMover.fwd.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>
#include <protocols/docking/util.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

// Core headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/excn/Exceptions.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/pointer/deep_copy.hh>

#include <string>

namespace protocols {
namespace glycan_docking {


/// @brief TODO
/// @author Morgan L. Nance (@mlnance)
class GlycanDockProtocol : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default GlycanDock constructor
	GlycanDockProtocol();

	/// @brief Copy constructor for deep copies
	GlycanDockProtocol ( GlycanDockProtocol const & object_to_copy );

	/// @brief Assignment operator
	//GlycanDockProtocol & operator=( GlycanDockProtocol const & object_to_copy );

	/// @brief Destructor
	// important for properly forward-declaring smart-pointer members
	~GlycanDockProtocol() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Apply the GlycanDockProtocol mover
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Generate string representation of GlycanDockProtocol
	void
	show( std::ostream & output ) const override;

	/// @brief PyRosetta-friendly version of show
	void
	show() const;


	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag
	// to use GlycanDockProtocol Mover in Rosetta Scripts
	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/////////////////////
	/// Setter Methods //
	/////////////////////

	/// @brief Set the ScoreFunction used for the GlycanDockProtocol
	/// GlycanDock was benchmarked using the ref2015 scoring function
	/// with sugar_bb = 0.5 and fa_intra_rep_nonprotein 0.55 (fa_rep weight)
	void
	set_scorefunction( core::scoring::ScoreFunctionOP const & in_sf ) { scorefxn_ = in_sf; }

	/// @brief Set the number of times to repeat the GlycanDockProtocol
	/// @brief starting from the decoy conformation after Stage 1 (if applied)
	/// @brief if the decoy fails the docking filter at the end of the protocol
	void
	set_n_repeats( core::Size const & n_repeats ) { n_repeats_ = n_repeats; }

	/// @brief Set whether to skip Stage 1 of GlycanDockProtocol and perform more
	/// @brief conservative glycosidic torsion angle moves during Stage 2
	void
	set_refine_only( bool const refine_only ) { refine_only_ = refine_only; }

	/// @brief Set whether to apply only the GlycanDock pre-packing protocol
	void
	set_prepack_only( bool const prepack_only ) { prepack_only_ = prepack_only; }

	/// @brief Set the chain IDs identifying the protein-glycoligand complex
	/// @brief protein is chain A, glycoligand is chain X --> partners = A_X
	void
	set_partners( std::string const & p ) { partners_ = p; }

	/// @brief Set whether the Stage 1 conformation initialization portion
	/// @brief of the GlycanDockProtocol is performed
	void
	set_stage1_rotate_glycan_about_com( bool const sg1 ) { stage1_rotate_glyc_about_com_ = sg1; }

	/// @brief Set the magnitude used when applying a random translational
	/// @brief perturbation to the initial glycoligand conformation in Stage 1
	void
	set_stage1_perturb_glycan_com_trans_mag( core::Real const trans_mag ) { stage1_trans_mag_ = trans_mag; }

	/// @brief Set the magnitude used when applying a random rotational
	/// @brief perturbation to the initial glycoligand conformation in Stage 1
	void
	set_stage1_perturb_glycan_com_rot_mag( core::Real const rot_mag ) { stage1_rot_mag_ = rot_mag; }

	/// @brief Set the magnitude used when applying a random perturbation to
	/// @brief each initial glycosidic torsion angle during Stage 1
	void
	set_stage1_torsion_uniform_pert_mag( core::Real const pert_mag ) { stage1_torsion_uniform_pert_mag_ = pert_mag; }

	/// @brief Set the number of rigid-body sampling and optimization
	/// @brief rounds to perform during each inner cycle of Stage 2
	void
	set_n_rigid_body_rounds( core::Size const rb_rounds ) { n_rb_rounds_ = rb_rounds; }

	/// @brief Set wehter to ocassionally rigidly slide the glycoligand
	/// @brief toward the center-of-mass of its proteein docking partner
	/// @brief during Stage 2 rigid-body sampling and optimization cycles
	void
	set_slide_glycan_into_contact( bool const slide ) { slide_glyc_into_contact_ = slide; }

	/// @brief Set the number of glycosidic torsion angle sampling and optimization
	/// @brief rounds to perform during each inner cycle of Stage 2
	void
	set_n_torsion_rounds( core::Size const tor_rounds ) { n_tor_rounds_ = tor_rounds; }

	/// @brief Set the magnitude used when applying a random translational
	/// @brief perturbation to the glycoligand conformation in Stage 2
	void
	set_stage2_trans_mag( core::Real const trans_mag ) { stage2_trans_mag_ = trans_mag; }

	/// @brief Set the magnitude used when applying a random rotational
	/// @brief perturbation to the glycoligand conformation in Stage 2
	void
	set_stage2_rot_mag( core::Real const rot_mag ) { stage2_rot_mag_ = rot_mag; }

	/// @brief Set whether to use a random, non-branch-point residue of the
	/// @brief glycoligand to set as the Jump residue for the docking FoldTree
	void
	set_rand_glycan_jump_res( bool const rand_glyc_jump ) { rand_glyc_jump_res_ = rand_glyc_jump; }

	/// @brief Set the temperature (kT) used when applying the Metropolis
	/// @brief criterion via the MonteCarlo object to accept or a reject a move
	void
	set_mc_kt( core::Real const mc_kt ) { mc_kt_ = mc_kt; }

	/// @brief Set the number of times to apply the more costly PackRotamersMover
	/// @brief during the inner cycles of Stage 2 sampling and optimization.
	/// @details Should be a factor of n_rb_rounds and n_tor_rounds. When the
	/// @details PackRotamersMover is not used, the EnergyCutRotamerTrialsMover is
	void
	set_full_packing_frequency( core::Size const full_pack_freq ) { full_pack_freq_ = full_pack_freq; }

	/// @brief Set the distance in Ang. defining the protein-glycoligand interface
	/// @brief when performing packing during Stage 2 sampling and optimization
	/// @note Uses less accurate neighbor atom distances to get interface
	void
	set_interface_packing_distance( core::Real const intf_pack_dist ) { intf_pack_dist_ = intf_pack_dist; }

	/// @brief Set whether to ramp the fa_atr and fa_rep score terms of the
	/// @brief ScoreFunction used during Stage 2 sampling and optimization
	void
	set_ramp_scorefxn( bool const ramp_sf ) { ramp_scorefxn_ = ramp_sf; }

	/// @brief Send the conformation of the decoy to PyMOL at specific,
	/// @brief hard-coded points during the GlycanDockProtocol
	/// @note TODO - currently not implemented
	void
	set_watch_in_pymol( bool const pymol ) { watch_in_pymol_ = pymol; }

	/////////////////////
	/// Getter Methods //
	/////////////////////

	/// @brief Get the ScoreFunction used for the GlycanDockProtocol
	core::scoring::ScoreFunctionOP const
	get_scorefunction() { return scorefxn_; }

	/// @brief Get the number of times to repeat the GlycanDockProtocol
	/// @brief starting from the decoy conformation after Stage 1 (if applied)
	/// @brief if the decoy fails the docking filter at the end of the protocol
	core::Size
	get_n_repeats() { return n_repeats_; }

	/// @brief Get whether to skip Stage 1 of GlycanDockProtocol and perform more
	/// @brief conservative glycosidic torsion angle moves during Stage 2
	bool
	get_refine_only() { return refine_only_; }

	/// @brief Get whether to apply only the GlycanDock pre-packing protocol
	bool
	get_prepack_only() { return prepack_only_; }

	/// @brief Get the chain IDs identifying the protein-glycoligand complex
	/// @brief protein is chain A, glycoligand is chain X --> partners = A_X
	std::string
	get_partners() const { return partners_; }

	/// @brief Get whether the Stage 1 conformation initialization portion
	/// @brief of the GlycanDockProtocol is performed
	bool
	get_stage1_rotate_glycan_about_com() { return stage1_rotate_glyc_about_com_; }

	/// @brief Get the magnitude used when applying a random translational
	/// @brief perturbation to the initial glycoligand conformation in Stage 1
	core::Real
	get_stage1_perturb_glycan_com_trans_mag() { return stage1_trans_mag_; }

	/// @brief Get the magnitude used when applying a random rotational
	/// @brief perturbation to the initial glycoligand conformation in Stage 1
	core::Real
	get_stage1_perturb_glycan_com_rot_mag() { return stage1_rot_mag_; }

	/// @brief Get the magnitude used when applying a random perturbation to
	/// @brief each initial glycosidic torsion angle during Stage 1
	core::Real
	get_stage1_torsion_uniform_pert_mag() { return stage1_torsion_uniform_pert_mag_; }

	/// @brief Get the number of rigid-body sampling and optimization
	/// @brief rounds to perform during each inner cycle of Stage 2
	core::Size
	get_n_rigid_body_rounds() { return n_rb_rounds_; }

	/// @brief Get wehter to ocassionally rigidly slide the glycoligand
	/// @brief toward the center-of-mass of its proteein docking partner
	/// @brief during Stage 2 rigid-body sampling and optimization cycles
	bool
	get_slide_glycan_into_contact() { return slide_glyc_into_contact_; }

	/// @brief Get the number of glycosidic torsion angle sampling and optimization
	/// @brief rounds to perform during each inner cycle of Stage 2
	core::Size
	get_n_torsion_rounds() { return n_tor_rounds_; }

	/// @brief Get the magnitude used when applying a random translational
	/// @brief perturbation to the glycoligand conformation in Stage 2
	core::Real
	get_stage2_trans_mag() { return stage2_trans_mag_; }

	/// @brief Get the magnitude used when applying a random rotational
	/// @brief perturbation to the glycoligand conformation in Stage 2
	core::Real
	get_stage2_rot_mag() { return stage2_rot_mag_; }

	/// @brief Get whether to use a random, non-branch-point residue of the
	/// @brief glycoligand to set as the Jump residue for the docking FoldTree
	bool
	get_rand_glycan_jump_res() { return rand_glyc_jump_res_; }

	/// @brief Get the temperature (kT) used when applying the Metropolis
	/// @brief criterion via the MonteCarlo object to accept or a reject a move
	core::Real
	get_mc_kt() { return mc_kt_; }

	/// @brief Get the number of times to apply the more costly PackRotamersMover
	/// @brief during the inner cycles of Stage 2 sampling and optimization.
	/// @details Should be a factor of n_rb_rounds and n_tor_rounds. When the
	/// @details PackRotamersMover is not used, the EnergyCutRotamerTrialsMover is
	core::Size
	get_full_packing_frequency() { return full_pack_freq_; }

	/// @brief Get the distance in Ang. defining the protein-glycoligand interface
	/// @brief when performing packing during Stage 2 sampling and optimization
	/// @note Uses less accurate neighbor atom distances to get interface
	core::Real
	get_interface_packing_distance() { return intf_pack_dist_; }

	/// @brief Get whether to ramp the fa_atr and fa_rep score terms of the
	/// @brief ScoreFunction used during Stage 2 sampling and optimization
	bool
	get_ramp_scorefxn() { return ramp_scorefxn_; }

	/// @brief Send the conformation of the decoy to PyMOL at specific,
	/// @brief hard-coded points during the GlycanDockProtocol
	/// @note TODO - currently not implemented
	bool
	get_watch_in_pymol() { return watch_in_pymol_; }


private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register options from command line upon construction
	void register_options();

	/// @brief Setup defaults not managed by flags at construction
	void set_defaults();

	/// @brief Initialize additional command-line options specified by user
	void set_cmd_line_defaults();

	/// @brief Set protocol values from command line upon construction
	void init_options();

	/// @brief Copy protocol data
	//void
	//copy_data( GlycanDockProtocol & object_to_copy_to,
	// GlycanDockProtocol const & object_to_copy_from );


	////////////////////////
	/// Protocol Methods ///
	////////////////////////

	/// @brief Reset counters used for diagnostic purposes during docking
	void clear_counters();

	/// @brief Prepare scoring function with necessary sugar score terms
	/// Also stores target weights for fa_atr and fa_rep for ramping
	void setup_scorefxn( core::pose::Pose & pose);

	/// @brief Setup function for preparing FoldTree for glycoligand docking
	/// By default, a random (non branch point) carbohydrate residue
	/// will be selected as the downstream Jump residue
	void setup_docking_foldtree
	( core::pose::Pose & pose,
		core::pose::Pose & ref_pose );

	/// @brief Set the weight for the specified ScoreType in the ScoreFunction
	/// Score ramping based on ramping factor and fraction of protocol completion
	///@note Assumes that only fa_atr and fa_rep are being ramped
	void ramp_score_weight
	( core::scoring::ScoreType const score_type,
		core::Real const target_weight,
		core::Size const current_cycle,
		core::Size const n_cycles );

	/// @brief Stage 2 docking refinement portion of the GlycanDockProtocol
	/// Performs inner cycles of n_rb_rounds_ of rigid-body perturbations
	/// and n_tor_rounds_ of glycosidic torsion angle perturbation
	/// Calling this function is considered one inner cycle of Stage 2
	/// Packing is performed after every perturbation, where the ecut_packer
	/// is called every time the full_packer is not, and the full_packer is
	/// is called after every full_pack_req_ moves
	/// Example: n_rb_rounds_ = n_tor_rounds_ = 8; full_pack_freq_ = 8
	/// Minimization is called after every other perturbation
	void perform_stage2_docking_and_optimization
	( core::pose::Pose & pose,
		protocols::moves::RandomMoverOP stage2_rb_randmover,
		protocols::moves::RandomMoverOP stage2_tor_randmover,
		minimization_packing::PackRotamersMoverOP full_packer,
		minimization_packing::EnergyCutRotamerTrialsMoverOP ecut_packer,
		minimization_packing::MinMoverOP min_mover,
		moves::MonteCarloOP mc );

	/// @brief Check decoy score to determine if it passed the docking filter
	/// Default filter is checking if the decoy has an interaction_energy < 0
	bool docking_filter( core::pose::Pose const & pose,
		core::Size const repeat );

	/// @brief Record a collection of decoy metrics
	void record_pose_metrics
	( core::pose::Pose & pose,
		std::string const & partners,
		utility::vector1< bool > const & glycolig_subset,
		core::pose::Pose const & ref_pose,
		utility::vector1< bool > const & ref_glycolig_subset );

	/// @brief Relay and record protocol-level GlycanDock results
	void record_protocol_info( core::pose::Pose & pose );

	/// @brief Prepack protocol for GlycanDock (based on FlexPepDock prepack)
	///         (1)  translate glycoligand away (by 1000 Ang)
	///         (2)  independently pack all side-chains
	///         (3)  translate glycoligand back to starting position
	void prepack_only
	( core::pose::Pose & pose,
		protocols::minimization_packing::PackRotamersMoverOP prepack_packer );


private: // data

	//////////////////////////
	//////// Scoring /////////
	//////////////////////////

	utility::pointer::DeepCopyOP< core::scoring::ScoreFunction > scorefxn_;

	core::Real mc_kt_ = 0.6;
	bool ramp_scorefxn_ = true; // cmd line
	core::Real target_atr_ = 1; // must be filled by starting fa_atr weight
	core::Real target_rep_ = 0.55; // must be filled by starting fa_rep weight
	core::Real starting_ramp_down_atr_factor_ = 3.25; // values by Krishna
	core::Real starting_ramp_up_rep_factor_ = 0.45455; // values by Krishna


	//////////////////////////
	//////// Movers //////////
	//////////////////////////

	// Packing
	core::Size full_pack_freq_ = 8; // cmd line
	core::Real const rt_energycut_ = 0.05; // from benchmark and FlexPepDock

	// Minimization
	bool lock_rings_ = true; // cmd line
	std::string min_type_ = "lbfgs_armijo_nonmonotone"; // cmd line; from benchmark
	core::Real min_tolerance_ = 0.01; // cmd line; from benchmark


	//////////////////////////
	//// Protocol Behavior ///
	//////////////////////////

	core::Size n_cycles_ = 10; // from benchmark, set by -run:ncycles
	core::Size n_repeats_ = 1; // cmd line; benchmark used = 3
	std::string partners_ = "_"; // req
	core::Real intf_pack_dist_ = 16; // cmd line

	bool rand_glyc_jump_res_ = true; // cmd line
	core::Size dock_jump_num_ = 1;
	protocols::docking::DockJumps movable_jump_;

	bool prepack_only_ = false; // cmd line
	bool refine_only_ = false; // cmd line
	bool watch_in_pymol_ = false; // cmd line


	//////////////////////////
	// Selectors and Subsets /
	//////////////////////////
	// We assume the input glycoligand is at least a disaccharide or larger
	// if input glycoligand is a monosaccharide, it has no torsions to sample
	// in that case, we would skip performing torsion sampling and instead
	// perform additional rigid-body sampling
	bool glycoligand_has_dihedrals_ = true;


	//////////////////////////
	//// Stage 1 Constants ///
	//////////////////////////

	core::Real stage1_trans_mag_ = 0.5; // cmd line
	core::Real stage1_rot_mag_ = 7.5; // cmd line
	bool stage1_rotate_glyc_about_com_ = false; // cmd line
	core::Real stage1_torsion_uniform_pert_mag_ = 12.5; // cmd line


	//////////////////////////
	//// Stage 2 Constants ///
	//////////////////////////

	core::Size n_rb_rounds_ = 8; // cmd line
	core::Real stage2_trans_mag_ = 0.5; // cmd line
	core::Real stage2_rot_mag_ = 7.5; // cmd line
	bool slide_glyc_into_contact_ = true; // cmd line
	core::Size n_tor_rounds_ = 8; // cmd line
	core::Size n_shear_moves_ = 1; // for ShearMover


	//////////////////////////
	//// General Constants ///
	//////////////////////////

	// counters for post-protocol analysis
	core::Size n_rb_moves_made_ = 0;
	core::Size n_rb_moves_accepted_ = 0;
	core::Size n_rb_cycles_performed_ = 0;
	core::Size n_tor_moves_made_ = 0;
	core::Size n_tor_moves_accepted_ = 0;
	core::Size n_tor_cycles_performed_ = 0;
	core::Size n_stage2_cycles_performed_ = 0; // debug; ensure it == n_cycles

};


std::ostream &
operator<<( std::ostream & os, GlycanDockProtocol const & mover );

} //glycan_docking
} //protocols

#endif //protocols_glycan_docking_GlycanDockProtocol_HH
