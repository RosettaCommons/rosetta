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

#ifndef INCLUDED_core_scoring_constraints_CoordinateConstraint_hh
#define INCLUDED_core_scoring_constraints_CoordinateConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Conformation.fwd.hh>


// C++ Headers
#include <cstdlib>
#include <iostream>
//#include <map>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


/// @details CoordinateConstraint compares the coordinates of a given atom (AtomID atom_) to a
/// fixed coordinate triplet (Vector xyz_target_).  Its other argument, fixed_atom_, is somewhat
/// nonobvious.  CoordinateConstraints are meant to be used with a Pose that has a nonmoving
/// virtual root residue.  An AtomID in this virtual root residue should be passed as fixed_atom_.
/// CoordinateConstraint does not use fixed_atom_, but the ScoreFunction code detects when
/// fixed_atom_ and atom_ move relative to one another, and trigger re-scoring at that time.  In
/// other words, CoordinateConstraints are really context-independent one body energies, but we
/// wish them to be evaluated as context-independent two-body energies.  (Ideally, ScoreFunction
/// would detect when atom_ moves relative to xyz_target_, but since ScoreFunction functions on
/// atoms and not floating coordinate triplets, this is a good workaround.) -- SML
class CoordinateConstraint : public Constraint {
public:


	CoordinateConstraint();

	/// @brief Constructor as a deep copy of the src %CoordinateConstraint
	CoordinateConstraint( CoordinateConstraint const & src );

	///c-tor
	CoordinateConstraint(
		AtomID const & a1,
		AtomID const & fixed_atom_in,
		Vector const & xyz_target_in,
		core::scoring::func::FuncOP func,
		ScoreType scotype = coordinate_constraint
	);

	~CoordinateConstraint();


	virtual std::string type() const;

	virtual ConstraintOP clone() const;

	virtual ConstraintOP clone( core::scoring::func::FuncOP ) const;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=NULL ) const;


	void show( std::ostream& out ) const;

	// @brief Reads the definition of a Constraint from the given std::istream,
	// using the given Pose, and the given func::FuncFactory. This method is intended
	// to be overridden by derived classes if they'd like to use the
	// ConstraintIO machinery.
	virtual void read_def( std::istream &, pose::Pose const &,func::FuncFactory const & );

	/// @brief possibility to do object comparison instead
	/// of pointer comparison
	virtual bool operator == ( Constraint const & other_cst) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	// @brief take coordinates, distances, angles, etc from given pose
	///
	virtual void steal_def( pose::Pose const& );

	Real
	non_virtual_score(
		Vector const & xyz
	) const;


	virtual
	void
	score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	// atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;


	Size
	natoms() const;

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;


	AtomID const &
	atom( Size const n ) const;

	// UNDEFIUNED, commenting out to fix PyRosetta build  void set_atom( Size const index, AtomID const atom_id );

	Real
	dist( pose::Pose const & pose ) const;


	virtual Size show_violations(
		std::ostream& out,
		pose::Pose const& pose,
		Size verbose_level,
		Real threshold = 1
	) const;

private:

	// functions
	Real
	func( Real const theta ) const;

	// deriv
	Real
	dfunc( Real const theta ) const;

private:
	// data
	AtomID atom_;
	AtomID fixed_atom_;
	Vector xyz_target_;
	core::scoring::func::FuncOP func_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_CoordinateConstraint )
#endif // SERIALIZATION


#endif
