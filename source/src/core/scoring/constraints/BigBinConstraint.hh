// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/BigBinConstraint.hh
/// @author James Thompson

#ifndef INCLUDED_core_scoring_constraints_BigBinConstraint_hh
#define INCLUDED_core_scoring_constraints_BigBinConstraint_hh

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/func/FuncFactory.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>


#include <core/id/AtomID.hh>

//Utility Headers

// C++ Headers
#include <iostream>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


/// constraint on dihedral angle formed by 4 points

class BigBinConstraint : public Constraint {
public:

	virtual std::string type() const;

	virtual ConstraintOP clone() const;

	virtual bool operator == ( Constraint const & other ) const;
	virtual bool same_type_as_me( Constraint const & other ) const;

	void score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const;

	// atom deriv
	void
	fill_f1_f2(
		AtomID const &,
		func::XYZ_Func const & xyz,
		Vector &,
		Vector &,
		EnergyMap const &
	) const;

	///c-tor
	BigBinConstraint(
		AtomID C0,
		AtomID N1,
		AtomID CA1,
		AtomID C1,
		AtomID N2,
		AtomID CA2,
		char bin,
		ScoreType scotype = dihedral_constraint
	);
	BigBinConstraint();
	BigBinConstraint( Size const res, char const bin, core::Real const sdev );


	Size
	natoms() const;

	AtomID const &
	atom( Size const n ) const;

	virtual void show( std::ostream & out ) const;

	void read_def( std::istream & in, pose::Pose const & pose,func::FuncFactory const & func_factory );

	char bin() const {
		return bin_;
	}

	Size res() const {
		return res_;
	}

	Real sdev() const {
		return sdev_;
	}
protected:
	BigBinConstraint( BigBinConstraint const & src );

private:
	// data
	AtomID C0_, N1_, CA1_, C1_, N2_, CA2_;
	Size res_;
	char bin_;
	core::Real sdev_;

	utility::vector1< ConstraintOP > my_csts_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class BigBinConstraint

typedef utility::pointer::shared_ptr< BigBinConstraint > BigBinConstraintOP;
typedef utility::pointer::shared_ptr< BigBinConstraint const > BigBinConstraintCOP;

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_BigBinConstraint )
#endif // SERIALIZATION


#endif
