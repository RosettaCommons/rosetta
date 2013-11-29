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

#ifndef INCLUDED_core_scoring_constraints_AtomPairConstraint_hh
#define INCLUDED_core_scoring_constraints_AtomPairConstraint_hh

#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

#include <utility/vector1.hh>



// C++ Headers
//#include <cstdlib>
//#include <iostream>
//#include <map>
//#include <utility>


namespace core {
namespace scoring {
namespace constraints {


///

class AtomPairConstraint : public Constraint {
public:

	// default c-tor
	AtomPairConstraint() : Constraint( atom_pair_constraint ) {}

	///c-tor
	AtomPairConstraint(
		AtomID const & a1,
		AtomID const & a2,
	 	FuncOP func,
		ScoreType scoretype = atom_pair_constraint
	):
		Constraint( scoretype ),
		atom1_(a1),
		atom2_(a2),
		func_( func )
	{}

	virtual ConstraintOP clone() const {
		return new AtomPairConstraint( atom1_, atom2_, func_, score_type() );
	}

	///
	virtual
	ConstraintOP clone( FuncOP func ) const {
		return new AtomPairConstraint( atom1_, atom2_, func, score_type() );
	}


	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	bool operator == ( Constraint const & other ) const;

	using Constraint::score;

	///
	Real
	score(
		Vector const & xyz1,
		Vector const & xyz2
		) const;

	///
	void
	score( XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	Real score( pose::Pose const& pose ) const {
		return func_->func( dist( pose ) );
	}

	// atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		XYZ_Func const & xyz,
		Vector & F1,
	 	Vector & F2,
		EnergyMap const & weights
	) const;

	std::string type() const {
		return "AtomPair";
	}

	///
	Size
	natoms() const
	{
		return 2;
	}

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;

	///
	AtomID const &
	atom( Size const n ) const
	{
		switch( n ) {
		case 1:
			return atom1_;
		case 2:
			return atom2_;
		default:
			utility_exit_with_message( "AtomPairConstraint::atom() bad argument" );
		}
		return atom1_;
	}

	void show( std::ostream& out ) const;
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	void read_def( std::istream& in, pose::Pose const& pose, FuncFactory const& func_factory );
	//	//@brief set constraint such that the pose doesn't violate it.
	//	virtual void steal( pose::Pose& );
	virtual
	Real dist( pose::Pose const& pose ) const;

	Real dist( core::conformation::Conformation const& conformation ) const;

	virtual
	Real dist( XYZ_Func const & xyz ) const;

	virtual Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const;

	virtual Func const& get_func() const {
		return *func_;
	}

	virtual
	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const;


private:
	// functions
	Real
	func( Real const theta ) const
	{
		return func_->func( theta );
	}

	// deriv
	Real
	dfunc( Real const theta ) const
	{
		return func_->dfunc( theta );
	}

protected:
	// data -- write accessed by NamedAtomPairConstraint
	mutable AtomID atom1_, atom2_;
	FuncOP func_;
}; // class AtomPairConstraint

} // constraints
} // scoring
} // core

#endif
