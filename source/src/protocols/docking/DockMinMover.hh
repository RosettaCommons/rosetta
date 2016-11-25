// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file docking_min_protocol
/// @brief
/// @author Robin A Thottungal (raugust1@jhu.edu)
#ifndef INCLUDED_protocols_docking_DockMinMover_hh
#define INCLUDED_protocols_docking_DockMinMover_hh

// Unit Headers
#include <protocols/docking/DockMinMover.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/TrialMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace docking {

//DockMinMover
class DockMinMover : public DockingHighRes {
public:
	/// @brief Default constructor
	DockMinMover();

	/// @brief Constructor with two arguments.  The first is the DockJumps, the second is a scorefunction to use for
	/// minimization.
	DockMinMover(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn
	);
	//JQX: constructor with mc_ object
	DockMinMover(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn,
		moves::MonteCarloOP mc
	);


	/// @brief Constructor with seven arguments.  The first is the DockJumps, the second is a movemap, the third is a
	/// scorefunction, the fourth is the minimization type (as a string), the fifth is the tolerance, the sixth is a boolean
	/// for the nb_list and the seventh is a MonteCarloOP
	DockMinMover(
		DockJumps const movable_jumps,
		core::kinematics::MoveMapOP movemap,
		core::scoring::ScoreFunctionOP scorefxn,
		std::string min_type,
		core::Real min_tolerance,
		bool nb_list,
		moves::MonteCarloOP mc
	);

	~DockMinMover() override;

	/// @brief Sets up the default values for the obejct including the movemap and minimization type.
	void set_default();

	/// @brief setters for member variables
	void set_min_type( std::string min_type ) { min_type_ = min_type; }
	void set_min_tolerance( core::Real min_tolerance ) { min_tolerance_ = min_tolerance; }

	void apply( core::pose::Pose & ) override;

	// XRW TEMP  std::string get_name() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	core::kinematics::MoveMapOP movemap_;
	std::string min_type_;
	core::Real min_tolerance_;
	bool nb_list_;
	moves::MonteCarloOP mc_;
	moves::TrialMoverOP minimize_trial_;
};

}
}
#endif
