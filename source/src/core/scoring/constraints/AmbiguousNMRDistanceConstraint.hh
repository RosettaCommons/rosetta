// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_core_scoring_constraints_AmbiguousNMRDistanceConstraint_hh
#define INCLUDED_core_scoring_constraints_AmbiguousNMRDistanceConstraint_hh

#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>

#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/chemical/AA.hh>

#include <utility/vector1.hh>


// C++ Headers
//#include <cstdlib>
//#include <iosfwd>
//#include <map>
//#include <utility>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


class AmbiguousNMRDistanceConstraint : public Constraint {
public:
	typedef utility::vector1< AtomID > Atoms;

public:
	/// @brief c-tor
	AmbiguousNMRDistanceConstraint(
		Atoms const & a1,
		Atoms const & a2,
		func::FuncOP func,
		ScoreType scoretype = atom_pair_constraint
	);

	AmbiguousNMRDistanceConstraint(
		id::NamedAtomID const & a1, //digests names like "QG1"
		id::NamedAtomID const & a2,
		core::pose::Pose const&,
		func::FuncOP func,
		ScoreType scoretype = atom_pair_constraint
	);

	AmbiguousNMRDistanceConstraint();

	ConstraintOP clone() const override;

	ConstraintOP clone( func::FuncOP func ) const override;

	bool operator == ( Constraint const & other ) const override;

	bool same_type_as_me( Constraint const & other ) const override;

	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. nullptr = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	ConstraintOP remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP map=nullptr ) const override;

	void read_def( std::istream& in, pose::Pose const& pose,func::FuncFactory const& func_factory ) override;
	// //@brief set constraint such that the pose doesn't violate it.
	// virtual void steal( pose::Pose& );

	///returns AtomPairConstraint or AmbigousNMRDistanceConstraint (e.g. for GLY HA1-HA2 ... )
	ConstraintOP map_to_CEN( pose::Pose const& src, pose::Pose const& centroid, core::Size& nr_mapped, std::string const& map_atom ) const;

	// Needed to get the base class overloads
	using Constraint::score;
	using Constraint::dist;

	void
	score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const override;

	core::Real
	inv_dist6( func::XYZ_Func const & xyz ) const;

	Real dist( func::XYZ_Func const & xyz ) const override;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const override;

	std::string type() const override {
		return "AmbiguousNMRDistance";
	}


	Size
	natoms() const override
	{
		return atoms1_.size() + atoms2_.size();
	}


	Size
	natoms( core::Size i ) const
	{
		return i==1 ? atoms1_.size() : atoms2_.size();
	}

	Size
	multiplicity() const {
		return atoms1_.size()*atoms2_.size();
	}

	/// @brief return residue number: i=1,2
	Size
	resid( core::Size i ) const {
		return i == 1 ? atoms1_[ 1 ].rsd() : atoms2_[ 1 ].rsd();
	}

	ConstraintOP
	remap_resid( core::id::SequenceMapping const &seqmap ) const override;


	AtomID const &
	atom( Size const n ) const override
	{
		if ( n > natoms() ) {
			utility_exit_with_message( "AmbiguousNMRDistanceConstraint::atom() bad argument" );
		}
		if ( n <= atoms1_.size() ) {
			return atoms1_[ n ];
		} else return atoms2_[ n - atoms1_.size() ];
	}

	void show( std::ostream& out ) const override;
	void show_def( std::ostream& out, pose::Pose const& pose ) const override;

	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const override;

	func::Func const& get_func() const override {
		return *func_;
	}

	core::Size effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const override;

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
	Atoms const & atoms1() const;
	Atoms const & atoms2() const;

	void atoms1( Atoms const & setting );
	void atoms2( Atoms const & setting );

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	AmbiguousNMRDistanceConstraint( AmbiguousNMRDistanceConstraint const & src );

private:
	Atoms atoms1_, atoms2_;
	func::FuncOP func_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class AmbiguousNMRDistanceConstraint

typedef utility::vector1< core::id::NamedAtomID > NamedAtoms;
void parse_NMR_name( std::string name, core::Size res, core::chemical::AA aa, NamedAtoms& );

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_AmbiguousNMRDistanceConstraint )
#endif // SERIALIZATION


#endif
