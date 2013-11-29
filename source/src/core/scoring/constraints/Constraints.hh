// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_constraints_Constraints_hh
#define INCLUDED_core_scoring_constraints_Constraints_hh

// Unit headers
#include <core/scoring/constraints/Constraints.fwd.hh>

// Package headers
#include <core/scoring/constraints/Constraint.fwd.hh>
#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif
#include <core/scoring/func/XYZ_Func.fwd.hh>

/// Project headers
#include <core/types.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>

// Utility Headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class Constraints : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Constraints();
	typedef id::AtomID AtomID;
	typedef conformation::Residue Residue;
	typedef conformation::Conformation Conformation;
	typedef ConstraintCOPs::const_iterator const_iterator;

public:

	Constraints();
	Constraints( Constraints const & );
	ConstraintsOP clone() const;
	Constraints const & operator = ( Constraints const & );

	// call the setup_for_derivatives for each constraint
 	void
 	setup_for_scoring( XYZ_Func const & xyz_func, ScoreFunction const &scfxn ) const;

	// call the setup_for_derivatives for each constraint
 	void
 	setup_for_derivatives( XYZ_Func const & xyz_func, ScoreFunction const &scfxn ) const;

	/// will fail if Residues dont contain all the necessary atoms
	void
	residue_pair_energy(
		Residue const & rsd1,
		Residue const & rsd2,
		EnergyMap const & weights,
 		EnergyMap & emap
	) const;

	/// will fail if Residue doesnt contain all the necessary atoms
	void
	intra_residue_energy(
		Residue const & rsd,
		EnergyMap const & weights,
		EnergyMap & emap
	) const;

	void
	conformation_energy(
		Conformation const & conformation,
		EnergyMap const & weights,
		EnergyMap & emap
	) const;

	/// @brief Evaluate derivatives giving the Constraint objects held within this object
	/// a single residue.  Warning: if this Constraints object contains Constraint objects
	/// that operate on other residues besides the one being given them, then this function
	/// will cause the program to exit.
	void
	eval_intrares_atom_derivative(
		id::AtomID const & atom_id,
		conformation::Residue const & residue,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	/// @brief Evaluate derivatives giving the Constraint objects held within this object
	/// two residues.  Warning: if this Constraints object contains Constraint objects
	/// that operate on other residues besides the two being given them, then this function
	/// will cause the program to exit.
	void
	eval_respair_atom_derivative(
		id::AtomID const & atom_id,
		conformation::Residue const & residue1,
		conformation::Residue const & residue2,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	/// @brief Evaluate derivatives giving the Constraint objects held within this object
	/// the entire Conformation (a whole structure, ws) with which to work.
 	void
 	eval_ws_atom_derivative(
 		AtomID const & atom_id,
 		Conformation const & conformation,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	///
	void
	add_constraint( ConstraintCOP cst );

	/// @brief Returns true iff the constraint was successfully found and removed.
	bool
	remove_constraint(
		ConstraintCOP cst,
		bool object_comparison
	);

	const_iterator begin() const;
	const_iterator end() const;

	void
	show( std::ostream& out );

	void
	show_definition( std::ostream& out, pose::Pose const& pose ) const;

	virtual Size
	show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, core::Real threshold = 1  );

	Size
	size() const;

	void
	clear();

	ConstraintCOPs const& constraints() const { return constraints_; }

private:
	void
	energy( XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const;

	// There's no implementation for this method in the .cc file ...
	//void
	//add_residue_pair_constraint( Size const pos1, Size const pos2, ConstraintCOP cst );

	void
	copy_from( Constraints const & );

private:

	ConstraintCOPs constraints_;

};


} // constraints
} // scoring
} // core

#endif
