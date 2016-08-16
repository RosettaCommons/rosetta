// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Roland A. Pache, PhD

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_FragmentPerturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_FragmentPerturber_HH

// Core headers
#include <core/fragment/FragSet.hh>

// Utility headers
#include <utility/vector1.hh>

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Pick phi/psi torsions from the given fragment library.
class FragmentPerturber : public Perturber {

public:

	/// @brief Constructor setting the fragment library.
	FragmentPerturber(utility::vector1< core::fragment::FragSetCOP > const & frag_libs);

	/// @copydoc Perturber::get_name
	std::string get_name() const { return "FragmentPerturber"; }

	/// @copydoc Perturber::get_subset
	void perturb_subset(
		Pose const & pose,
		IndexList const & residues,
		ClosureProblemOP problem);

private:

	/// @brief Fragment library.
	utility::vector1< core::fragment::FragSetCOP > frag_libs_;

};

}
}
}

#endif

