// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBP_FinalFullatomRefinementMoveGenerator.hh
/// @brief A class to generate ParsedProtocols for the final full-atom refinement step.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_helical_bundle_predict_HBP_FinalFullatomRefinementMoveGenerator_hh
#define INCLUDED_protocols_helical_bundle_predict_HBP_FinalFullatomRefinementMoveGenerator_hh

#include <protocols/helical_bundle_predict/HBP_MoveGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_FinalFullatomRefinementMoveGenerator.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief A class to generate ParsedProtocols for the final full-atom refinement step.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HBP_FinalFullatomRefinementMoveGenerator : public HBP_MoveGenerator {

public:

	/// @brief Default constructor is explicitly deleted.
	HBP_FinalFullatomRefinementMoveGenerator() = delete;

	/// @brief Options constructor.
	/// @details The provided scorefunction is cloned.
	HBP_FinalFullatomRefinementMoveGenerator( core::Size const fast_relax_rounds, bool const find_disulfides, core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Copy constructor.
	HBP_FinalFullatomRefinementMoveGenerator(HBP_FinalFullatomRefinementMoveGenerator const & src);

	/// @brief Destructor.
	~HBP_FinalFullatomRefinementMoveGenerator() override;

	/// @brief Clone this object: that is, make a copy and return an owning pointer to the copy.
	HBP_MoveGeneratorOP clone() const override;

public: //Functions

	/// @brief Given the current step index, the total number of steps in the trajectory, and the pose for analysis,
	/// construct a ParsedProtocol of things to do to this pose for this move in the Monte Carlo trajectory.
	/// @details Pure virtual.  Must be implemented by derived classes.
	protocols::rosetta_scripts::ParsedProtocolOP
	generate_monte_carlo_move(
		core::Size const current_step,
		core::Size const num_steps,
		core::pose::Pose const &pose
	) const override;

private: //Data

	/// @brief Number of FastRelax rounds to perform.
	core::Size fast_relax_rounds_;

	/// @brief Should we try disulfide permutations first?
	bool find_disulfides_;

	/// @brief The fullatom scorefunction used by FastRelax.
	core::scoring::ScoreFunctionOP sfxn_;

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HBP_FinalFullatomRefinementMoveGenerator_hh





