// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.hh
/// @brief In stepwise design, vary desired loop lengths by updating FullModelParameters
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_stepwise_monte_carlo_mover_VaryLoopLengthMover_hh
#define INCLUDED_protocols_stepwise_monte_carlo_mover_VaryLoopLengthMover_hh

// Unit headers
#include <protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

///@brief In stepwise design, vary desired loop lengths by updating FullModelParameters
class VaryLoopLengthMover : public protocols::moves::Mover {

public:

	VaryLoopLengthMover();

	// copy constructor
	VaryLoopLengthMover( VaryLoopLengthMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~VaryLoopLengthMover();

	void
	apply( core::pose::Pose & pose,  StepWiseMove const & swa_move );

	using protocols::moves::Mover::apply;
	virtual void
	apply( core::pose::Pose & pose );

public:

	std::string
	get_name() const;

private:

	void
	update_full_model_parameters( core::pose::Pose & pose,
																core::pose::full_model_info::FullModelParametersCOP full_model_parameters_new ) const;


};



} //protocols
} //stepwise
} //monte_carlo
} //mover


#endif //protocols/stepwise/monte_carlo/mover_VaryLoopLengthMover_hh







