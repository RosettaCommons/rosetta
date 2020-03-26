// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.hh
/// @brief A class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory, based on the
/// current state of the pose.  This version uses helix-coil transition theory to nucleate and extend helices.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_helical_bundle_predict_HBP_HelixCoilMoveGenerator_hh
#define INCLUDED_protocols_helical_bundle_predict_HBP_HelixCoilMoveGenerator_hh

#include <protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.fwd.hh>
#include <protocols/helical_bundle_predict/HBP_MoveGenerator.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>
#include <protocols/helical_bundle_predict/HBPHelixAssignments.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief A class for a module to generate ParsedProtocols for the next move in a Monte Carlo trajectory, based on the
/// current state of the pose.  This version uses helix-coil transition theory to nucleate and extend helices.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBP_HelixCoilMoveGenerator : public HBP_MoveGenerator {

public:

	/// @brief Constructor.
	HBP_HelixCoilMoveGenerator();

	/// @brief Copy constructor.
	HBP_HelixCoilMoveGenerator(HBP_HelixCoilMoveGenerator const & src);

	/// @brief Destructor.
	~HBP_HelixCoilMoveGenerator() override;

	/// @brief Clone this object: that is, make a copy and return an owning pointer to the copy.
	HBP_MoveGeneratorOP clone() const override;

public: //Functions

	/// @brief Given the current step index, the total number of steps in the trajectory, and the pose for analysis,
	/// construct a ParsedProtocol of things to do to this pose for this move in the Monte Carlo trajectory.
	/// @details This override uses helix-coil transition theory to nucleate, extend, or remove helix segments.  Non-helix
	/// segments are simply subjected to small moves.
	/// @returns A ParsedProtocolOP with the new protocol for success, or nullptr for failure.
	protocols::rosetta_scripts::ParsedProtocolOP
	generate_monte_carlo_move(
		core::Size const current_step,
		core::Size const num_steps,
		core::pose::Pose const &pose
	) const override;

	/// @brief When a move is accepted, set the current helix assignments to the candidate helix assignments.
	void mark_move_accepted() const override;

	/// @brief Given a helix assignment file's contents, set up the helix assignments.
	/// @details Involves no read from disk.
	void set_up_user_helix_assignments( std::string const &file_contents );

private: //Functions

	/// @brief Update the helix assignments to nucleate helices.
	void add_helix_nucleation_moves( core::pose::Pose const &pose ) const;

	/// @brief Update the helix assignments to elongate helices.
	void add_helix_elongation_moves( core::pose::Pose const &pose ) const;

	/// @brief Update the helix assignments to shrink helices.
	void add_helix_retraction_moves( core::pose::Pose const &pose ) const;

	/// @brief Add global parameter perturbations.
	void add_global_parameter_perturbations( core::pose::Pose const &pose ) const;

	/// @brief Add local parameter perturbations.
	void add_local_parameter_perturbations( core::pose::Pose const &pose ) const;

	/// @brief Add movers to the move list to update the helices based on the helix assignments.
	/// @details Returns true for FAILURE (due to inability to generate helix from parameters), false for success.
	bool add_helix_update_moves( core::pose::Pose const & pose, protocols::rosetta_scripts::ParsedProtocol &protocol ) const;

	/// @brief Add the SmallMover moves to the protocol, at non-helix positions.
	void add_nonhelix_small_moves( core::pose::Pose const &pose, protocols::rosetta_scripts::ParsedProtocol &protocol ) const;


private: //Data

	/// @brief Stores information about what parts of a pose might be helical, as provided
	/// by the user.
	HBPHelixAssignments user_helix_assignments_;

	/// @brief Stores information about what parts of a pose are currently believed to be helical, at a
	/// given step in a trajectory.
	mutable HBPHelixAssignments current_helix_assignments_;

	/// @brief Tentatively-considered helix assignments.
	mutable HBPHelixAssignments candidate_helix_assignments_;

	/// @brief A RamaPrePro scorefunction for general use of SmallMover.
	core::scoring::ScoreFunctionOP ramaprepro_sfxn_;

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBP_HelixCoilMoveGenerator_hh





