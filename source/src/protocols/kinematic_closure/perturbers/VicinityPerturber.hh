// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_VicinityPerturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_VicinityPerturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/VicinityPerturber.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Limit sampling to phi/psi pairs within some distance of a target 
/// loop.
class VicinityPerturber : public Perturber {

public:

	/// @brief Constructor which takes a target.  Sampling will be limited to 
	/// phi/psi pairs within some delta of this target.
	VicinityPerturber(Pose const & target);

	/// @copydoc Perturber::get_name
	string get_name() const { return "VicinityPerturber"; }

	/// @copydoc Perturber::perturb_subset
	void perturb_subset(
			Pose const & pose,
			IndexList const & residues,
			ClosureProblemOP problem);

	/// @copydoc Perturber::perturb_subset_with_balance
	void perturb_subset_with_balance(
			Pose const & pose,
			IndexList const & residues,
			ClosureProblemOP problem);

private:
	Pose const target_;
	Real const spread_;

};

}
}
}

#endif

