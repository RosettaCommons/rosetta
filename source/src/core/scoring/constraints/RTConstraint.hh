// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/RTConstraint.hh
/// @brief Constraint definition for RT of stub pairs
/// @author atom-moyer (apmoyer@uw.edu)

#ifndef INCLUDED_core_scoring_constraints_RTConstraint_hh
#define INCLUDED_core_scoring_constraints_RTConstraint_hh

// package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>

#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>

#include <core/kinematics/RT.hh>
#include <core/conformation/Conformation.fwd.hh>

// C++ Headers
#include <cstdlib>
#include <iostream>
//#include <map>

#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief Actually a *restraint*, like a virtual rubber band between a pair of atoms.
///
/// @details All Constraints are expected to be immutable once created,
/// meaning their internal data (state) should not change over their lifetime.
/// This allows Constraints to be shared between copies of Poses (e.g. in Monte Carlo),
/// and is important for both speed (with thousands of constraints) and correctness.
///
/// To "change" a constraint, remove the old one and add a new and different one.
/// The steal() methods have been removed because it is
/// incompatible with the idea of immutable constraints.

class RTConstraint;
typedef utility::pointer::shared_ptr< RTConstraint > RTConstraintOP;
typedef utility::pointer::shared_ptr< RTConstraint const > RTConstraintCOP;

class RTConstraint : public Constraint {
public:

	RTConstraint();

	// @brief Constructor
	RTConstraint(
		id::StubID const & stub1,
		id::StubID const & stub2,
		kinematics::RT const & rt_target,
		func::FuncOP func,
		ScoreType const & score_type = atom_pair_constraint
	);

	// @brief Copy Constructor
	RTConstraint( RTConstraint const & src );

	/// @brief destructor
	~RTConstraint();

	/// @brief Copies the data from this %Constraint into a new object and returns
	/// an OP to the new object. Intended to be implemented by derived classes and
	/// used by pose.add_constraint.  This function must return a *deep copy* of
	/// itself -- meaning that if this %Constraint holds pointers to other %Constraints
	/// that it must invoke clone on those %Constraints as well.  If the %Constraint
	/// holds a FuncOP, then the Func should also be cloned.
	ConstraintOP
	clone() const override;

	/// @brief Clone the constraint, but where a new Func object is to be used instead.
	ConstraintOP
	clone( func::FuncOP ) const override;

	/// @brief calculates the RT from stub1->stub2 at this instant
	kinematics::RT
	calculate_rt( func::XYZ_Func const & ) const;

	/// @brief return the raw "distance" before that distance is handed to the FUNC object
	Real
	dist( func::XYZ_Func const & xyz ) const override;

	/// @brief Calculates a score for this constraint using XYZ_Func, and puts
	/// the UNWEIGHTED score into emap. Although the current set of weights
	/// currently is provided, Constraint objects should put unweighted scores
	/// into emap because the ScoreFunction will do the weighting itself.
	void
	score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const override;

	std::string
	type() const override;

	Size
	natoms() const override;

	AtomID const &
	atom( Size const n ) const override;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const override;

	bool
	operator == ( Constraint const & rhs ) const override;

	bool
	same_type_as_me( Constraint const & other ) const override;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	ConstraintOP
	remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=NULL ) const override;

	void
	show( std::ostream& out ) const override;

	func::Func
	const & get_func() const override;

private: // Data
	id::StubID stub1_;
	id::StubID stub2_;
	kinematics::RT rt_target_;
	func::FuncOP func_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class RTConstraint

typedef utility::pointer::shared_ptr< RTConstraint > RTConstraintOP;
typedef utility::pointer::shared_ptr< RTConstraint const > RTConstraintCOP;


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_RTConstraint )
#endif // SERIALIZATION


#endif
