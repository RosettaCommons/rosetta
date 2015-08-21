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
/// @details Solving a kinematic closure problem has two steps.  The first step
/// is to pick a new torsion angles, bond angles, and/or bond lengths in the
/// region being sampled.  The second step is to analytically set six torsions
/// such that the backbone stays closed.  Perturber subclasses are responsible
/// for managing the first step.
///
/// Every perturber subclass must reimplement perturb_subset().  The arguments
/// to the method are a Pose, a list of residues, and a ClosureProblem.  The
/// pose is const and should just be used to look up residue types and other
/// relevant contextual information.  The residue list specifies the residues
/// that should be perturbed.  In other words, residues not in this list should
/// not be changed.  This allows composite perturbers to specify different
/// perturber algorithms for different parts on the loop.  For example, for
/// antibody modeling you might use a custom perturber for loop regions with
/// well-known motifs and the standard RamaPerturber everywhere else.  Finally,
/// the problem is used to actually set the backbone DOFS via methods like
/// perturb_phi() and perturb_psi().
///
/// If your perturber algorithm can be made to obey detailed balance, you
/// should also reimplement perturb_subset_with_balance().  This method is
/// called by BalancedKicMover and it works just like perturb_subset().  If
/// your algorithm doesn't obey detailed balance, or if you don't know what
/// detailed balance is, then don't worry about it.
///
/// @note Currently perturb_phi() and perturb_psi() live in ClosureProblem.
/// This makes sense, because the problem owns the matrices that these methods
/// are perturbing.  So you can think of these methods as part of an interface
/// that allows the closure problem to be defined.  However, this approach will
/// not work well with non-canonical backbones.  An alternative is to expose a
/// much more general interface to the problem (i.e. get_torsion_angles() and
/// set_torsion_angles()).  Then subclasses of Perturber that are only meant to
/// work with one sort of backbone can be created, and these subclasses can
/// define methods like perturb_phi() and perturb_psi() which use the more
/// general ClosureProblem interface.  In any case, you should be aware that
/// the perturber interface may change in the near future.
///
/// It would have been nice if any regular mover could be used to perturb the
/// backbone, instead of being limited to single-purpose Perturber subclasses.
/// I spent a lot of time trying to get this to work, but it was very slow.  My
/// belief is that too many AtomTree coordinates were being updated too often.

class Perturber
	: public utility::pointer::ReferenceCount, protected boost::noncopyable {

public:

	/// @brief Return the name of this perturber.
	virtual std::string get_name() const = 0;

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
