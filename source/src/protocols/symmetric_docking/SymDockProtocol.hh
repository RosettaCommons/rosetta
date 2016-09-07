// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/symmetric_docking/SymDockProtocol.hh
///
/// @brief    Symmetric Protein-Protein Docking Protocol
/// @details  Dock symmetrical assemblies together - works with the symmetry
///           framework in Rosetta 3. Also includes some options for working with
///           membranes
///           Last Modified: 10/25/14
///
/// @author Ingemar Andre
/// @author Rebecca Alford (adding membranes & comments)

#ifndef INCLUDED_protocols_symmetric_docking_SymDockProtocol_hh
#define INCLUDED_protocols_symmetric_docking_SymDockProtocol_hh

// Unit Headers
#include <protocols/symmetric_docking/SymDockProtocol.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/symmetric_docking/SymDockingLowRes.fwd.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.fwd.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

// Package Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace symmetric_docking {

/// @brief Define Main funciton for symmetric docking protocol
void SymDock_main();

/// @brief Main Symmetric Docking protocol class
class SymDockProtocol : public moves::Mover
{
public:

	typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;

public:


	/// @brief Default Constructor for Symmetric Docking Protocol
	SymDockProtocol();

	/// @brief Custom Constructor - Symmetric docking protocol
	/// @details Custom protocol - specify if input pose is fullatom, local refinement
	/// and use graphics viewer
	SymDockProtocol(
		bool const fullatom,
		bool const local_refine,
		bool const view=false
	);

	/// @brief Custom Constructor - Symmetric Docking Protocol
	/// @details Custom protocol - specify if input pose is fullatom, local refinement
	/// and use graphics viewer. Also specify custom scoring functions for low and high resolution
	/// stages
	SymDockProtocol(
		bool const fullatom,
		bool const local_refine,
		bool const view,
		core::scoring::ScoreFunctionOP docking_score_low,
		core::scoring::ScoreFunctionOP docking_score_high
	);

	/// @brief Destructor
	~SymDockProtocol() override;

	/// @brief Register Options from the Commandline
	void register_options();

	/// @brief Setup options based on constructors
	void set_default();


	/// Setters for Protocol Options
	void set_dock_rtmin( bool dock_rtmin_in );

	void set_sc_min( bool sc_min_in );
	void set_max_repeats( Size const max_repeats_in );
	void set_dock_ppk( bool dock_ppk_in );

	void set_fullatom( bool const fullatom_in );

	void set_local_refine( bool const local_refine_in );

	void set_view( bool view_in );

	void set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_score_low_in );

	void set_highres_scorefxn( core::scoring::ScoreFunctionOP docking_score_high_in );

	void set_highres_scorefxn(
		core::scoring::ScoreFunctionOP docking_score_high_in,
		core::scoring::ScoreFunctionOP docking_score_pack_in
	);

	bool docking_lowres_filter( core::pose::Pose & pose );
	bool docking_highres_filter( core::pose::Pose & pose );

	/// @brief Calculate Interaction energy between partners
	core::Real
	calc_interaction_energy( core::pose::Pose & pose );

	/// @brief Calculate RMSD from partner to symmetric starting structure
	core::Real
	calc_rms( core::pose::Pose & pose );

	/// @brief Recover sidechains from the native pose
	void recover_sidechains( core::pose::Pose & pose, const core::pose::Pose & native_pose );

	/// @brief Set Task Factory for Packing
	void task_factory( core::pack::task::TaskFactoryOP task_factory );

	/// @brief Turn on Design of partner 2 during docking - not thoroughly tested
	/// @note @ralford I think this is just copied over from regular docking protocol?
	void design( bool const des );
	bool design() const;

	/// @brief Skip population of the score map
	void hurry( bool const hurry );

	/// @brief Methods for task factory?
	core::pack::task::TaskFactoryOP task_factory() const;
	core::pack::task::TaskFactoryOP & task_factory();

	/// @brief Only score the pose using docking config
	void score_only( core::pose::Pose & pose );

	/// RosettaScripts Methods //////
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/// Mover Methods /////////
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;

private:

	/// @brief Classic MCM Protocol for symmetric docking protocol
	void
	classic_mcm_protocol(
		core::pose::Pose & pose,
		core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo,
		core::Size num_cycles,
		core::Size repack_every_Nth
	) const;


	/// @brief Setup Docking Monte Carlo mover
	protocols::moves::MoverOP
	make_dockmcm_mover(
		core::pose::Pose const & pose,
		protocols::moves::MoverOP repack_mover,
		protocols::moves::MoverOP rigbod_mover,
		core::kinematics::MoveMapOP movemap, //< would be COP but MinMover wants OP
		core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo
	) const;

private:

	/// @brief Data about the protocol
	bool fullatom_;
	bool local_refine_;
	bool rtmin_;
	bool sc_min_;
	Size max_repeats_;
	bool dock_ppk_;

	/// @brief Jumps for symmetric rigid body transformations
	utility::vector1<core::Size> movable_jumps_;

	/// @brief Modified Docking foldtree for symmetry
	bool autofoldtree_;

	/// @brief Initialize the graphics viewer (opengl)
	bool view_;

	/// @brief Allow design of partner 2 during docking
	bool design_;

	/// @brief Filters for structures during protocols
	bool passed_lowres_filter_;
	bool passed_highres_filter_;

	/// @brief Sip population of the score map
	bool hurry_; // skip populating the score map

	/// @brief For scorefile output
	std::map < std::string, core::Real > score_map_;

	/// @brief Scorefunctions used during the protocol
	core::scoring::ScoreFunctionOP docking_score_low_;
	core::scoring::ScoreFunctionOP docking_score_high_;
	core::scoring::ScoreFunctionOP docking_score_high_min_;
	core::scoring::ScoreFunctionOP docking_score_pack_;

	/// @brief Monte carlo object
	moves::MonteCarloOP mc_;

	/// @brief Child protocols - low resolution and high resolution docking & for packing
	protocols::symmetric_docking::SymDockingLowResOP docking_low_;
	protocols::symmetric_docking::SymDockingHiResOP docking_high_;

	/// @brief Used to restrict packer task during docking
	core::pack::task::TaskFactoryOP init_task_factory_;

}; // SymDockProtocol

} // symmetric_docking
} // protocols

#endif //INCLUDED_protocols_symmetric_docking_SymDockProtocol_HH
