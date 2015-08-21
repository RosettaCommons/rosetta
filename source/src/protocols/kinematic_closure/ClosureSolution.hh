// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ClosureSolution.hh
/// @brief  Header file for ClosureSolution.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_kinematic_closure_ClosureSolution_HH
#define INCLUDED_protocols_kinematic_closure_ClosureSolution_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {

using utility::pointer::ReferenceCount;
using boost::noncopyable;
using core::pose::Pose;

/// @brief Represent a single solution to a kinematic closure problem.
///
/// @details The ClosureSolution class represents the solutions returned by
/// ClosureProblem.solve().  The most important methods of this class are
/// apply() and apply_if_reasonable().  The former unconditionally applies the
/// solution to the given pose, while the latter does so only if the solution
/// passes a rama and bump check.

class ClosureSolution : public ReferenceCount, private noncopyable {

	friend class ClosureProblem;

	// Constructors {{{1
private:

	/// @brief Constructor used internally to build a solution from the
	/// internal degrees of freedom returned by the closure algorithm.
	ClosureSolution(
		ClosureProblem const * problem,
		Size const solution_index,
		ParameterList const & torsion_angles,
		ParameterList const & bond_angles,
		ParameterList const & bond_lengths);

	// Interesting functions {{{1
public:

	/// @brief Apply this solution to the given pose.
	void apply(Pose & pose) const;

	/// @brief If this solution passes rama and bump checks, apply it to the
	/// given pose.  Return whether or not the filters were passed.
	bool apply_if_reasonable(
		Pose & pose,
		bool rama_on=true,
		bool bump_on=true,
		bool be_lenient=false) const;

	/// @brief Return a unique number identifying this solution.
	Size get_index() const;

	/// @brief Return the Jacobian for this solution.
	Real get_jacobian() const;

	/// @ brief Return a distance metric indicating how similar this solution
	/// is to the given problem.
	Real get_distance(ClosureProblem const * problem) const;

	// Private Helpers {{{1
private:

	/// @brief Check for unlikely pivot torsions in this solution.
	bool check_rama(Pose const & pose, Real const temperature) const;

	/// @brief Check for overlapping backbone atoms in this solution.
	bool check_overlap(Pose const & pose, Real const scale_factor) const;

	// Data Members {{{1
private:

	ClosureProblem const * problem_;
	Size const index_;

	ParameterList const bond_lengths_;
	ParameterList const bond_angles_;
	ParameterList const torsion_angles_;

	mutable Real jacobian_;

	// }}}1

};

} // namespace kinematic_closure
} // namespace protocols

#endif
