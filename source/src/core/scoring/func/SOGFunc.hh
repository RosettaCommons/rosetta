// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SOGFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_func_SOGFunc_hh
#define INCLUDED_core_scoring_func_SOGFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>

#include <core/scoring/func/SOGFunc.fwd.hh>
#include <core/scoring/func/SOGFunc_Impl.hh>

// C++ Headers

// AUTO-REMOVED #include <ostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace func {

/// @brief Derived class of class Func representing a SOG distribution with a user-specified
/// mean and standard deviation.
class SOGFunc : public Func {
public:

	SOGFunc() {}

	SOGFunc(
		const utility::vector1< core::Real >& means,
		const utility::vector1< core::Real >& sdevs,
		const utility::vector1< core::Real >& weights
	);

	SOGFunc( core::Real mean, core::Real sdev );

	/// @brief returns a clone of this SOGFunc
	FuncOP clone() const { return new SOGFunc( *this ); }

	/// @brief Returns the value of this SOGFunc evaluated at distance x.
	Real func( Real const x ) const;

	/// @brief Returns the value of the first derivative of this SOGFunc at distance x.
	Real dfunc( Real const x ) const;

	void check_bounds( Real const x, Real const val ) const;

	/// @brief show the definition of this SOGFunc to the specified output stream.
	virtual void show_definition( std::ostream &out ) const;

	/// @brief Calls show( out ) on this SOGFunc.
	friend std::ostream& operator<<(std::ostream& out, const SOGFunc& f ) {
		f.show( out );
		return out;
	} // operator<<


	/// @brief Initializes this SOGFunc from the given istream.
	/// @detailed The parameters are:
	/*!
	 * Initializes this SOGFunc from the given istream. An example
	 * of the type of string from which the istream should be constructed is:
	 * "SOGFUNC 2 19.396 7.643 0.4 17.312 2.4 0.6". The interpretation is to
	 * create initialize this SOGFunc object with the following parameters:
	 * - one Gaussian function with mean 19.396 and sd 7.643, weighted with a
	 * weight of 0.4
	 * - another Gaussian function with mean of 17.312 and sd 2.4, with a weight
	 * of 0.6.
	 * Weights need not add up to 1, but many times they will.
	 */
	void read_data( std::istream & in );

private:
	void clear_();
	core::Real get_alt_score_( core::Real const x ) const;

	SOGFunc_Impl member_;
};



} // constraints
} // scoring
} // core

#endif
