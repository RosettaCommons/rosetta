// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_BondAnglePerturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_BondAnglePerturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/BondAnglePerturber.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Resample each backbone bond angle from a gaussian distribution.
/// @details The parameters for the underlying gaussian distribution were taken
/// from a fit of the bond angle distribution observed in the Top8000 database
/// in 2013.  Also note that bond angles (3 atoms) are different than torsion
/// angles (4 atoms).

class BondAnglePerturber : public Perturber {

public:

	/// @copydoc Perturber::get_name()
	std::string get_name() const { return "BondAnglePerturber"; }

	/// @copydoc Perturber::perturb_subset()
	void perturb_subset(
		Pose const & pose,
		IndexList const & residues,
		ClosureProblemOP problem);

};

}
}
}

#endif
