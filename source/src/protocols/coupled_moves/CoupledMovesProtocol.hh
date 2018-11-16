// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/coupled_moves/CoupledMovesProtocol.hh
/// @brief Mover that implements the CoupledMovesProtocol
/// @author Noah <nollikai@gmail.com>, refactored into header by Steven Lewis smlewi@gmail.com
/// @author Anum Glasgow
/// @author Amanda Loshbaugh

#ifndef INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh
#define INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh

#include <protocols/coupled_moves/CoupledMovesProtocol.fwd.hh>
#include <protocols/simple_moves/CoupledMover.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/minimization_packing/BoltzmannRotamerMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinPackMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <ctime>
#include <iostream>
#include <string>
#include <sstream>

//fragment
#include <core/fragment/FragSet.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <protocols/loops/loops_main.hh>

//KIC
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>

namespace protocols {
namespace coupled_moves {

/// @namespace protocols::coupled_moves
///
/// @brief Protocol that combines backbone moves and sidechain moves in a single Monte Carlo step.
///
/// @details CoupledMoves protocol makes use of the simple mover CoupledMover. The CoupledMover apply function makes a single coupled move as follows:
/// (1) Backbone move -- As of March 2018, there are three options: backrub, kinematic closure with walking perturber, or kinematic closure with fragment perturber.
/// (2) Side chain move -- The three residues at the center of the perturbed backbone are sampled by BoltzmannRotamerMover.
/// (3) Repack -- 10A shell around the residues from step 2 is repacked.
///
/// CoupledMoves will also move a ligand.
///
/// Assumptions:
/// - If there is a ligand, it should be the last residue in the pdb.

class CoupledMovesProtocol : public protocols::moves::Mover {
public:
	CoupledMovesProtocol();
	CoupledMovesProtocol(CoupledMovesProtocol const & cmp);

	virtual void apply( core::pose::Pose& pose );

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

	///@brief adds terms for the backrub step to the scorefunction.  Used by ctor; also useful when using set_score_fxn
	void configure_score_fxn();
	core::Real compute_ligand_score_bonus( core::pose::PoseOP pose );
	void exclude_nonclashing_positions( core::pose::PoseOP pose );
	//void setup_move_positions( core::pack::task::PackerTaskOP task, core::pose::PoseOP pose );
	void setup_move_and_design_positions( core::pack::task::PackerTaskCOP task, core::pose::PoseCOP pose );
	protocols::simple_moves::CoupledMoverOP setup_coupled_mover( core::pack::task::PackerTaskOP task, core::pose::PoseOP pose_copy, protocols::simple_moves::CoupledMoverOP coupled_mover );
	protocols::simple_moves::CoupledMoverOP setup_coupled_mover_for_ligand( core::pose::PoseOP pose_copy, core::pack::task::PackerTaskOP task, protocols::simple_moves::CoupledMoverOP coupled_mover );
	void setup_move_type( core::Real const move_prob );
	void setup_resnum( core::Real const move_prob, protocols::simple_moves::CoupledMoverOP coupled_mover, core::pose::PoseOP pose_copy );
	void write_info_to_log( core::pose::PoseOP pose_copy );
	///@brief parse_my_tag
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) /*override*/;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	// setters
	///@brief set the score function.  NOTICE: you will want to also use the companion configure_score_fxn afterwards to add backrub terms to the SF, if you did not do so yourself!
	void set_score_fxn( core::scoring::ScoreFunctionOP const sf ) { score_fxn_ = sf; }
	///@brief set the task factory.  Note that the ctor sets InitializeFromCommandline plus either ReadResfile or RestrictToRepacking.
	void set_main_task_factory( core::pack::task::TaskFactoryOP const tf ) { main_task_factory_ = tf; }
	void set_ntrials( core::Size ntrials ) { ntrials_ = ntrials; }
	void set_loop_size( core::Size loop_size ) { loop_size_ = loop_size; }
	void set_perturber( kinematic_closure::perturbers::PerturberOP perturber ) { perturber_ = perturber; }
	void set_backbone_mover( std::string const & backbone_mover ) { backbone_mover_ = backbone_mover; }
	void set_repack_neighborhood ( bool repack_neighborhood ) { repack_neighborhood_ = repack_neighborhood; }
	// ligand setters
	void set_ligand_mode( bool const ligand_mode ) { ligand_mode_ = ligand_mode; }
	void set_number_ligands( core::Size const number_ligands ) { number_ligands_ = number_ligands; }
	void set_ligand_weight( core::Real const ligand_weight ) { ligand_weight_ = ligand_weight; }
	void set_ligand_prob( core::Real const ligand_prob ) { ligand_prob_ = ligand_prob; }

	// getters
	std::string get_name() const;
	static std::string mover_name();
	core::Size get_ntrials() { return ntrials_; }
	std::string get_backbone_mover() const { return backbone_mover_; }
	const bool & get_repack_neighborhood() const { return repack_neighborhood_; }
	core::Size get_loop_size() const { return loop_size_; }
	kinematic_closure::perturbers::PerturberOP get_perturber() const { return perturber_; }
	std::string get_move_type() const { return move_type_; }
	std::map<std::string,core::Real> get_unique_sequences() const { return unique_sequences_; }
	utility::vector1<core::Size> get_move_positions() const { return move_positions_; }
	utility::vector1<core::Size> get_design_positions() const { return design_positions_; }
	// ligand getters
	bool get_ligand_mode() const { return ligand_mode_; }
	core::Size get_number_ligands() const { return number_ligands_; }
	core::Real get_ligand_weight() const { return ligand_weight_; }
	core::Real get_ligand_prob() const { return ligand_prob_; }
	utility::vector1<core::Size> get_ligand_resnums() const { return ligand_resnums_; }

private:
	///@brief sets up internal options based on options system
	void initialize_from_options();
	///@brief sets up internal options for ligand mode based on the options system
	void initialize_ligand_from_options();

	core::scoring::ScoreFunctionOP score_fxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;
	core::Size ntrials_ = 1000; //default copied from options system
	std::string move_type_ = "RESIDUE"; //default to the safer case
	utility::vector1<core::Size> move_positions_;
	utility::vector1<core::Size> design_positions_;
	std::map<std::string,core::Real> unique_sequences_;
	std::map<std::string,core::scoring::EnergyMap> unique_scores_;
	/// @brief parameter for KIC with walking perturber
	//core::Real walking_perturber_magnitude_;
	/// @brief loop size for KIC
	core::Size loop_size_ = 4;
	/// @brief mover used for side-chain moves
	//protocols::minimization_packing::BoltzmannRotamerMoverOP boltzmann_rotamer_mover_;
	/// @brief perturber used by kic
	kinematic_closure::perturbers::PerturberOP perturber_;
	/// @brief type of backbone move, e.g. backrub or kic
	std::string backbone_mover_ = "backrub";
	/// @brief Default false for legacy behavior (Ollikainen 2015). After the backbone move and rotamer move, repack sidechains within 5A of the design residue.
	bool repack_neighborhood_ = false;

	//Ligand options.  These were once handled only by command line flag.
	/// @brief if true, model protein ligand interaction
	bool ligand_mode_ = false;
	/// @brief number of ligands in the pose
	core::Size number_ligands_ = 1; //I think this is a silly default but it maintains existing behavior; the value is irrelevant if ligand_mode is not active
	/// @brief weight for residue - ligand interactions
	core::Real ligand_weight_ = 1.0;
	/// @brief probability of making a ligand move
	core::Real ligand_prob_ = 0.1;
	utility::vector1<core::Size> ligand_resnums_;

}; //CoupledMovesProtocol

} //coupled_moves
} //protocols

#endif //INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_HH
