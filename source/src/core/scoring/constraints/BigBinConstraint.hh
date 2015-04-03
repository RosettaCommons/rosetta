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


namespace core {
namespace scoring {
namespace constraints {


/// constraint on dihedral angle formed by 4 points

class BigBinConstraint : public Constraint {

public:

	virtual std::string type() const {
		return "BigBin";
	}

	virtual ConstraintOP clone() const {
		return ConstraintOP( new BigBinConstraint( res_, bin_, sdev_ ) );
	}

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
	):
		Constraint( scotype ),
		C0_( C0 ),
		N1_( N1 ),
		CA1_( CA1 ),
		C1_( C1 ),
		N2_( N2 ),
		CA2_( CA2 ),
		bin_(bin)
	{}

  BigBinConstraint() :
		Constraint( dihedral_constraint ),
		res_( 0 ),
		bin_( 'A' ),
		sdev_( 0.5 )
  {}

	BigBinConstraint( Size const res, char const bin, core::Real const sdev ) :
		Constraint( dihedral_constraint ),
		res_( res ),
		bin_( bin ),
		sdev_( sdev )
	{}


	Size
	natoms() const {
		return 6;
	}


	AtomID const &
	atom( Size const n ) const
	{
		switch( n ) {
			case 1:
				return C0_;
			case 2:
				return N1_;
			case 3:
				return CA1_;
			case 4:
				return C1_;
			case 5:
				return N2_;
			case 6:
				return CA2_;
			default:
				utility_exit_with_message( "BigBinConstraint::atom() bad argument" );
		}
		return C0_;
	}

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

private:
	// data
	AtomID C0_, N1_, CA1_, C1_, N2_, CA2_;
  Size res_;
  char bin_;
	core::Real sdev_;

	utility::vector1< ConstraintOP > my_csts_;
}; // class BigBinConstraint

} // constraints
} // scoring
} // core

#endif
