// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/ConstantFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author John Karanicolas

#ifndef INCLUDED_core_scoring_constraints_ConstantFunc_hh
#define INCLUDED_core_scoring_constraints_ConstantFunc_hh

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

// C++ Headers


namespace core {
namespace scoring {
namespace constraints {

/// @brief Derived class of class Func representing a Constant distribution with a user-specified
/// mean and standard deviation.
class ConstantFunc : public Func {
public:

	// Constuctor for ConstantFunc. Arguments to the constructor are:
	// return_val: (constant) return value
	ConstantFunc ( Real const return_val ) :
		return_val_( return_val )
	{}

	/// @brief returns a clone of this ConstantFunc
	FuncOP clone() const { return new ConstantFunc( *this ); }

	/// @brief Returns the value of this ConstantFunc evaluated at distance x.
	Real func( Real const ) const;

	/// @brief Returns the value of the first derivative of this ConstantFunc at distance x.
	Real dfunc( Real const ) const;

	/// @brief show the definition of this ConstantFunc to the specified output stream.
	virtual void show_definition( std::ostream & out ) const;

	/// @brief Calls show( out ) on this ConstantFunc.
	friend std::ostream & operator<<(std::ostream & out, const ConstantFunc & f ) {
		f.show( out );
		return out;
	} // operator<<

	void read_data( std::istream & in );

private:

	Real return_val_;
};

} // constraints
} // scoring
} // core

#endif
