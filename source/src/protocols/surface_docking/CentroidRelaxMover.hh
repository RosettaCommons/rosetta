// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CentroidRelaxMover.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (rathottungal@gmail.com)
/// @author Michael Pacella (mpacella88@gmail.com)

#ifndef INCLUDED_protocols_surface_docking_CentroidRelaxMover_hh
#define INCLUDED_protocols_surface_docking_CentroidRelaxMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/surface_docking/CentroidRelaxMover.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

// Utility Headers

// C++ Headers
#include <string>

namespace protocols {
namespace surface_docking {

class CentroidRelaxMover : public moves::Mover {

public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief Default constructor
	CentroidRelaxMover();

	/// @brief Copy constructor
	CentroidRelaxMover(CentroidRelaxMover const & src);

	//Destructor
	~CentroidRelaxMover();

	// Standard Rosetta methods ////////////////////////////////////////////////
	//General Methods
	/// @brief Register options with the options system
	// Undefined, commenting out to fix PyRosetta build static void register_options();

	/// @brief Generate string representation of CentroidRelaxMover
	void show(std::ostream & output=std::cout) const;

	//Insertion operator (overloaded so that CentroidRelaxMover can be "printed" in PyRosetta.
	// Undefined, commenting out to fix PyRosetta build friend std::ostream & operator<<(std::ostream & output, CentroidRelaxMover const & object_to_output);

	// Mover Methods
	/// @brief Return the name of the Mover
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Apply the corresponding move to <input_pose>
	virtual void apply( core::pose::Pose & pose );

	//Accessors/Mutators
	void set_nmoves(const core::Size setting);

	void set_inner_loop_cycles(const core::Size setting);

	void set_outer_loop_cycles(const core::Size setting);


private:

	// Private methods /////////////////////////////////////////////////////////
	void init();
	void copy_data(CentroidRelaxMover object_to_copy_to, CentroidRelaxMover object_to_copy_from);
	void setup_movers(const core::pose::Pose & pose);
	void setup_defaults();
	// Private data ////////////////////////////////////////////////////////////
	core::Size nmoves_; // number of positions at which to make moves


	//MC Settings
	core::Real kT_;

	// for scoring
	core::scoring::ScoreFunctionOP score_low_res_;

	// for movers

	std::string small_min_type_;
	std::string shear_min_type_;

	moves::MonteCarloOP  monte_carlo_;
	moves::TrialMoverOP small_trial_min_mover_;
	moves::TrialMoverOP shear_trial_min_mover_;

	core::Size inner_loop_cycles_;
	core::Size outer_loop_cycles_;

};


} // surfaceDockingProtocols
} // protocols

#endif //INCLUDED_protocols_surface_docking_CentroidRelaxMover_hh
