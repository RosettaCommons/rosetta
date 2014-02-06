// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_Perturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_Perturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

/// @brief Base class for all of the perturber algorithms.
/// @details Subclasses of Perturber are meant to be created by the programmer 
/// but called internally by samplers::KicMover.  If you want to implement a 
/// new perturber, all you have to do is reimplement perturb_subset().  The 
/// perturb() method simply calls perturb_subset() with a complete list of 
/// residues.  If you're algorithm has an implementation which obeys detailed 
/// balance, you can also reimplement perturb_subset_with_balance().

class Perturber
	: public utility::pointer::ReferenceCount, protected boost::noncopyable {

public:

	/// @brief Return the name of this perturber.
	virtual string get_name() const = 0;

	/// @brief Perturb all of the non-pivot residues.
	void perturb(
			Pose const & pose,
			ClosureProblemOP problem);

	/// @brief Perturb all of the non-pivot residues such that detailed balance 
	/// is obeyed.  
	void perturb_with_balance(
			Pose const & pose,
			ClosureProblemOP problem);

	/// @brief Perturb the given residues.
	virtual void perturb_subset(
			Pose const & pose,
			IndexList const & residues,
			ClosureProblemOP problem) = 0;

	/// @brief Perturb the given residues such that detailed balance is obeyed.  
	virtual void perturb_subset_with_balance(
			Pose const & pose,
			IndexList const & residues,
			ClosureProblemOP problem);

};

}
}
}

#endif

