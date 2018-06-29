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

#ifndef INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh
#define INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh

#include <protocols/coupled_moves/CoupledMovesProtocol.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace coupled_moves {

class CoupledMovesProtocol : public protocols::moves::Mover {
public:
	CoupledMovesProtocol();
	CoupledMovesProtocol(CoupledMovesProtocol const & cmp);

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

	core::Real compute_ligand_score_bonus(
		core::pose::PoseOP pose,
		utility::vector1<core::Size> ligand_resnums);

	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) /*override*/;

	static std::string mover_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	///@brief set the score function.  NOTICE: you will want to also use the companion configure_score_fxn afterwards to add backrub terms to the SF, if you did not do so yourself!
	void set_score_fxn( core::scoring::ScoreFunctionOP const sf ) { score_fxn_ = sf; }

	///@brief set the task factory.  Note that the ctor sets InitializeFromCommandline plus either ReadResfile or RestrictToRepacking.
	void set_main_task_factory( core::pack::task::TaskFactoryOP const tf ) { main_task_factory_ = tf; }

	///@brief adds terms for the backrub step to the scorefunction.  Used by ctor; also useful when using set_score_fxn
	void configure_score_fxn();

	bool get_ligand_mode() const { return ligand_mode_; }
	core::Size get_number_ligands() const { return number_ligands_; }
	core::Real get_ligand_weight() const { return ligand_weight_; }
	core::Real get_ligand_prob() const { return ligand_prob_; }

	void set_ligand_mode( bool const ligand_mode ) { ligand_mode_ = ligand_mode; }
	void set_number_ligands( core::Size const number_ligands ) { number_ligands_ = number_ligands; }
	void set_ligand_weight( core::Real const ligand_weight ) { ligand_weight_ = ligand_weight; }
	void set_ligand_prob( core::Real const ligand_prob ) { ligand_prob_ = ligand_prob; }

private:
	///@brief sets up internal options based on options system
	void initialize_from_options();

	///@brief sets up internal options for ligand mode based on the options system
	void initialize_ligand_from_options();

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;

	//Ligand options.  These were once handled only by command line flag.
	/// @brief if true, model protein ligand interaction
	bool ligand_mode_ = false;
	/// @brief number of ligands in the pose
	core::Size number_ligands_ = 1; //I think this is a silly default but it maintains existing behavior; the value is irrelevant if ligand_mode is not active
	/// @brief weight for residue - ligand interactions
	core::Real ligand_weight_ = 1.0;
	/// @brief probability of making a ligand move
	core::Real ligand_prob_ = 0.1;

}; //CoupledMovesProtocol

} //coupled_moves
} //protocols

#endif //INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_HH
