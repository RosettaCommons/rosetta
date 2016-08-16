// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FullatomRelaxMover.hh
/// @brief
/// @author Robin A Thottungal  (rathottungal@gmail.com)
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_FullatomRelaxMover_hh
#define INCLUDED_protocols_surface_docking_FullatomRelaxMover_hh

// Unit Headers
#include <protocols/surface_docking/FullatomRelaxMover.fwd.hh>
// Package headers
#include <protocols/surface_docking/SurfaceOrientMover.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
// Project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>
#include <protocols/docking/DockingInitialPerturbation.fwd.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace surface_docking {

class FullatomRelaxMover : public moves::Mover {

public:

	FullatomRelaxMover();

	FullatomRelaxMover(FullatomRelaxMover const & src);

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	//destructor
	~FullatomRelaxMover();

	void apply( core::pose::Pose & pose);

	virtual std::string get_name() const;

	void set_nmoves(const core::Size);

	std::string get_sol_secondary_struct();
	std::string get_ads_secondary_struct();

	void set_surface_contact_mover( protocols::docking::FaDockingSlideIntoContactOP surface_contact_mover );

	void set_surface_orient_mover( SurfaceOrientMoverOP surface_orient );

	void set_surface_parameters(protocols::surface_docking::SurfaceParametersOP surface_parameters);

private:

	void copy_data(FullatomRelaxMover object_to_copy_to, FullatomRelaxMover object_to_copy_from);

	void inner_loop_refinement( core::pose::Pose & pose );

	void outer_loop_refinement_solution( core::pose::Pose & pose );

	void outer_loop_refinement_adsorbed( core::pose::Pose & pose );

	void reorient_and_slide_into_surface( core::pose::Pose & pose );

	void dock_mcm_on_surface( core::pose::Pose & pose );

	void output_solution_state( core::pose::Pose & pose );

	void refinement_cycle( core::pose::Pose & pose );

	void setup_defaults();

	void setup_movers( const core::pose::Pose & pose );

	void set_smallmovesize(Size scale);

	void set_ljrepulsion_weight(core::Real weight_scale);

	void set_ecounter(core::Size ecount);

	void calc_secondary_struct(core::pose::Pose & pose);

	void set_secondary_struct(core::pose::Pose & pose);

	void reposition_above_surface(core::pose::Pose & pose);

	//members for smallTrialMove
	core::Real kT_;
	core::Size nmoves_;
	core::Size encounter_;
	core::Size encounter_cycle_;

	std::map< char, core::Real > angle_max_;
	// for scoring
	core::scoring::ScoreFunctionOP score_high_res_;
	// for movers

	std::string small_min_type_;
	std::string shear_min_type_;

	core::kinematics::MoveMapOP move_map_;

	moves::MonteCarloOP  monte_carlo_;

	simple_moves::SmallMoverOP small_mover_;
	simple_moves::MinMoverOP small_min_mover_;
	moves::SequenceMoverOP small_sequence_mover_;
	moves::TrialMoverOP small_trial_min_mover_;

	simple_moves::ShearMoverOP shear_mover_;
	simple_moves::MinMoverOP shear_min_mover_;
	moves::SequenceMoverOP shear_sequence_mover_;
	moves::TrialMoverOP shear_trial_min_mover_;

	//protocols::docking::DockingHighResOP dockingHigh_res_;
	std::string sol_sec_struct_;
	std::string ads_sec_struct_;
	std::string sec_struct_;

	protocols::docking::FaDockingSlideIntoContactOP surface_contact_mover_;
	protocols::surface_docking::SurfaceOrientMoverOP surface_orient_;
	protocols::surface_docking::SurfaceParametersOP surface_parameters_;
	docking::DockMCMProtocolOP dock_mcm_;
	core::Size outer_loop_cycles_;
	core::Size inner_loop_cycles_;
};


} // surface_docking
} // protocols

#endif
