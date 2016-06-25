// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/FromScratchMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_FromScratchMover_HH
#define INCLUDED_protocols_stepwise_monte_carlo_FromScratchMover_HH

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

class FromScratchMover: public protocols::moves::Mover {

public:

	//constructor
	FromScratchMover();

	//destructor
	~FromScratchMover();

public:

	using moves::Mover::apply;

	void
	apply( core::pose::Pose & pose,
		utility::vector1<Size> const & residues_to_instantiate_in_full_model_numbering ) const;

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler );

private:

	void
	update_full_model_info_and_switch_focus_to_new_pose( pose::Pose & pose, pose::Pose & new_pose, utility::vector1< Size > const & resnum ) const;

	void
	sample_by_swa( pose::Pose & pose, Size const sample_res ) const;

private:

	protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
