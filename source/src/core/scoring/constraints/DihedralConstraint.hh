// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/DihedralConstraint.hh
/// @brief

#ifndef INCLUDED_core_scoring_constraints_DihedralConstraint_hh
#define INCLUDED_core_scoring_constraints_DihedralConstraint_hh

#include <core/scoring/constraints/DihedralConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>

#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


/// constraint on dihedral angle formed by 4 points

class DihedralConstraint : public Constraint {
public:

	using Constraint::score; // for virtual Real score( Pose const & );

	Real
	score(
		Vector const & xyz1,
		Vector const & xyz2,
		Vector const & xyz3,
		Vector const & xyz4
	) const;

	virtual std::string type() const {
		return "Dihedral";
	}

	/// @breif Perform a deep clone by cloning the func_ pointer held by this %DihedralConstraint
	virtual ConstraintOP clone() const;

	Size show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold = 1 ) const;

	void read_def(
		std::istream & in,
		pose::Pose const & pose,
		func::FuncFactory const & func_factory
	);

	/// @brief possibility to compare constraint according to data
	/// and not just pointers
	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	void
	score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	Real
	score( conformation::Conformation const & conformation  ) const;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const & atom,
		func::XYZ_Func const & xyz,
		Vector & F1,
		Vector & F2,
		EnergyMap const & weights
	) const;

	///c-tor
	DihedralConstraint(
		AtomID const & a1,
		AtomID const & a2,
		AtomID const & a3,
		AtomID const & a4,
		func::FuncOP func,
		ScoreType scotype = dihedral_constraint
	):
		Constraint( scotype ),
		atom1_(a1),
		atom2_(a2),
		atom3_(a3),
		atom4_(a4),
		func_( func )
	{}

	Size
	natoms() const
	{
		return 4;
	}

	virtual
	ConstraintOP
	remap_resid( core::id::SequenceMapping const & seqmap ) const;

	AtomID const &
	atom( Size const n ) const;

	virtual void show( std::ostream & out ) const;

	void show_def( std::ostream& out, pose::Pose const& pose ) const;

	virtual func::Func const& get_func() const {
		return *func_;
	}

	// Coppied Remapped_clone from AtomPairConstraint
	/// @brief Copies the data from this Constraint into a new object and returns an OP
	/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
	/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
	/// to the new object. Intended to be implemented by derived classes.
	virtual ConstraintOP remapped_clone(
		pose::Pose const & src,
		pose::Pose const & dest,
		id::SequenceMappingCOP map = NULL
	) const;

private:
	Real
	func( Real const theta ) const;

	Real
	dfunc( Real const theta ) const;

	// this guy doesn't need member data --> static
	static
	void
	helper(
		Vector const & M,
		Vector const & v,
		Vector const & w,
		Vector & F1,
		Vector & F2
	);

	static
	void
	p1_cosine_deriv(
		Vector const & p1,
		Vector const & p2,
		Vector const & p3,
		Vector const & p4,
		Real & x,
		Vector & F1,
		Vector & F2
	);

	static
	void
	p2_cosine_deriv(
		Vector const & p1,
		Vector const & p2,
		Vector const & p3,
		Vector const & p4,
		Real & x,
		Vector & F1,
		Vector & F2
	);

protected:
	/// @brief Explicit copy constructor so that derived classes will recieve a deep copy
	/// of the Func this class contains.
	DihedralConstraint( DihedralConstraint const & src );


private:
	// data
	AtomID atom1_, atom2_, atom3_, atom4_;
	func::FuncOP func_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DihedralConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_DihedralConstraint )
#endif // SERIALIZATION


#endif
