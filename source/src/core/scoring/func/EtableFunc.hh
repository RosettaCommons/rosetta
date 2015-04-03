// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/EtableFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_func_EtableFunc_hh
#define INCLUDED_core_scoring_func_EtableFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>


// C++ Headers

#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace func {

/// @brief Class for representing arbitrarily defined functions.
class EtableFunc : public Func {
public:

	/*!
	 * Constuctor for EtableFunc. Arguments to the constructor are:
	 * - mean: parameter representing the mean of this function.
	 * - sd: parameter representing the standard deviation of this function.
	 */

	EtableFunc (
		Real const min_,
		Real const max_,
		Real const stepsize_
	) :
		min_     ( min_ ),
		max_     ( max_ ),
		stepsize_( stepsize_ )
	{}

	/// @brief returns a clone of this EtableFunc
	FuncOP clone() const { return FuncOP( new EtableFunc( *this ) ); }

	/// @brief Returns the value of this EtableFunc evaluated at distance x.
	Real func( Real const x ) const;

	/// @brief Returns the value of the first derivative of this EtableFunc at distance x.
	Real dfunc( Real const x ) const;

	/// @brief show the definitio of this EtableFunc to the specified output stream.
	virtual void show_definition( std::ostream &out ) const;

	/// @brief Calls show( out ) on this EtableFunc.
	friend std::ostream& operator<<(std::ostream& out, const EtableFunc& f ) {
		f.show( out );
		return out;
	} // operator<<


	void read_data( std::istream& in );

private:
	utility::vector1< core::Real > func_;
	core::Real min_;
	core::Real max_;
	core::Real stepsize_;
};


} // constraints
} // scoring
} // core

#endif
