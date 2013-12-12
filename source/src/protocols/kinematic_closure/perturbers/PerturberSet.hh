// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_PerturberSet_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_PerturberSet_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Maintain a set of perturbers to be executed together.
class PerturberSet : public Perturber {

public:

	/// @brief Default constructor
	PerturberSet();

	/// @brief Destructor
	~PerturberSet();

	/// @copydoc Perturber::get_name
	string get_name() const { return "PerturberSet"; }

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

	/// @brief Add a new perturber to the set.  Any default perturbers are 
	/// removed.
	void add(PerturberOP perturber);

	/// @brief Remove all perturbers from this set.
	void clear();

	/// @brief Indicate that the current set of perturbers is meant as some sort 
	/// of default, and should be cleared if a new perturber is manually added.
	void mark_as_default();

private:
	utility::vector1<PerturberOP> perturbers_;
	bool is_default_;

};

}
}
}

#endif

