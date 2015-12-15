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
#include <core/scoring/func/Func.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>

#include <utility/vector1.hh>


// C++ Headers
//#include <cstdlib>
//#include <iostream>
//#include <map>
//#include <utility>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class AtomPairConstraint : public Constraint {
public:

	// default c-tor
	AtomPairConstraint();

	///c-tor
	AtomPairConstraint(
		AtomID const & a1,
		AtomID const & a2,
		core::scoring::func::FuncOP func,
		ScoreType scoretype = atom_pair_constraint
	);

	/// @brief Create a deep copy of this %AtomPairConstraint, cloning its Func
	virtual ConstraintOP clone() const;

	/// @brief Create a deep copy of this %AtomPairConstraint except that the copy should
	/// use the input Func instead of its existing one.
	virtual
	ConstraintOP clone( func::FuncOP func ) const;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

	/// @brief Compare a) the class types (w/ same_type_as_me), b) the atoms being constrained,
	/// c) the score_type being used, and d) the Func objects (the FuncOPs may point at different
	/// objects, but as long as those objects are equal, that counts)
	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	using Constraint::score;


	Real
	score(
		Vector const & xyz1,
		Vector const & xyz2
	) const;


	void
	score( core::scoring::func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	Real score( pose::Pose const& pose ) const;

	// atom deriv
	virtual
	void
	fill_f1_f2(
		AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	std::string type() const {
		return "AtomPair";
	}


	Size
	natoms() const
	{
		return 2;
	}

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const;


	AtomID const &
	atom( Size const n ) const
	{
		switch( n ) {
		case 1 :
			return atom1_;
		case 2 :
			return atom2_;
		default :
			utility_exit_with_message( "AtomPairConstraint::atom() bad argument" );
		}
		return atom1_;
	}

	AtomID const & atom1() const { return atom1_; }
	AtomID const & atom2() const { return atom2_; }

	void show( std::ostream& out ) const;
	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	void read_def( std::istream& in, pose::Pose const& pose, func::FuncFactory const & func_factory );
	// //@brief set constraint such that the pose doesn't violate it.
	// virtual void steal( pose::Pose& );
	virtual
	Real dist( pose::Pose const& pose ) const;

	Real dist( core::conformation::Conformation const& conformation ) const;

	virtual
	Real dist( core::scoring::func::XYZ_Func const & xyz ) const;

	virtual Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const;

	virtual func::Func const & get_func() const;

	virtual
	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const;

private:
	// functions
	Real
	func( Real const theta ) const;

	// deriv
	Real
	dfunc( Real const theta ) const;

protected:
	/// @brief Setter for the derived class
	void set_func( func::FuncOP setting );

	/// @brief Setter for the mutable atom1_ data member to be used by
	/// derived classes.  Remember: data must never be protected, only
	/// private.  Instead, protected mutators can be added if derived
	/// classes ought to have the ability to change the base class data
	/// that other classes/functions should not.
	void atom1( AtomID newid ) const;

	/// @brief Setter for the mutable atom2_ data member to be used by
	/// derived classes.  Remember: data must never be protected, only
	/// private.  Instead, protected mutators can be added if derived
	/// classes ought to have the ability to change the base class data
	/// that other classes/functions should not.
	void atom2( AtomID newid ) const;

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	AtomPairConstraint( AtomPairConstraint const & src );

private:
	// no data
	mutable AtomID atom1_, atom2_;
	core::scoring::func::FuncOP func_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class AtomPairConstraint

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_AtomPairConstraint )
#endif // SERIALIZATION


#endif
