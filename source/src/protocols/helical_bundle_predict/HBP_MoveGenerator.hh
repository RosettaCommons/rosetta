// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_MoveGenerator.hh
/// @brief A base class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory, based on the current state of the pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_helical_bundle_predict_HBP_MoveGenerator_hh
#define INCLUDED_protocols_helical_bundle_predict_HBP_MoveGenerator_hh

#include <protocols/helical_bundle_predict/HBP_MoveGenerator.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief A base class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory,
/// based on the current state of the pose.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBP_MoveGenerator : public utility::VirtualBase {

public:

	/// @brief Constructor.
	HBP_MoveGenerator();

	/// @brief Copy constructor.
	HBP_MoveGenerator(HBP_MoveGenerator const & src);

	/// @brief Destructor.
	~HBP_MoveGenerator() override;

	/// @brief Clone this object: that is, make a copy and return an owning pointer to the copy.
	virtual HBP_MoveGeneratorOP clone() const = 0;

public: //Functions

	/// @brief Given the current step index, the total number of steps in the trajectory, and the pose for analysis,
	/// construct a ParsedProtocol of things to do to this pose for this move in the Monte Carlo trajectory.
	/// @details Pure virtual.  Must be implemented by derived classes.
	virtual
	protocols::rosetta_scripts::ParsedProtocolOP
	generate_monte_carlo_move(
		core::Size const current_step,
		core::Size const num_steps,
		core::pose::Pose const &pose
	) const = 0;

	/// @brief Does nothing in base class.  Can be overridden in derived classes to do something
	/// more useful when a move is accepted.
	virtual void mark_move_accepted() const;

	/// @brief Does nothing in base class.  Can be overridden in derived classes to do something
	/// more useful when a move is rejected.
	virtual void mark_move_rejected() const;

	/// @brief Set the current round.
	void set_current_round( core::Size const setting );

	/// @brief Set the total nunmber of rounds.
	void set_max_rounds( core::Size const setting );

	/// @brief Get the current round.
	inline core::Size current_round() const { return current_round_; }

	/// @brief Get the maximum number of rounds.
	inline core::Size max_rounds() const { return max_rounds_; }

private: //Data

	/// @brief This trajectory is which round in the series of trajectories?
	core::Size current_round_;

	/// @brief How many rounds (i.e. trajectories) are there ultimately?
	core::Size max_rounds_;

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBP_MoveGenerator_hh





