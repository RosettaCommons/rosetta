// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SOGFunc_Impl.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_constraints_SOGFunc_Impl_hh
#define INCLUDED_core_scoring_constraints_SOGFunc_Impl_hh

#include <core/scoring/func/SOGFunc_Impl.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

// C++ Headers
// AUTO-REMOVED #include <ostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

/// @brief Derived class of class Func representing a SOG distribution with a user-specified
/// mean and standard deviation.
class SOGFunc_Impl : public utility::pointer::ReferenceCount {
public:

	/*!
	 * Constuctor for SOGFunc_Impl. Arguments to the constructor are:
	 * - mean: parameter representing the mean of this function.
	 * - sd: parameter representing the standard deviation of this function.
	 */

	SOGFunc_Impl(
		const utility::vector1< core::Real >& means,
		const utility::vector1< core::Real >& sdevs,
		const utility::vector1< core::Real >& weights
	) :
		means_  ( means   ),
		sdevs_  ( sdevs   ),
		weights_( weights )
	{
		upper_bound(10);
	}

	SOGFunc_Impl() {
		means_.push_back  ( 0.0 );
		sdevs_.push_back  ( 1.0 );
		weights_.push_back( 1.0 );
		upper_bound(10);
	}
	virtual ~SOGFunc_Impl();

	/// @brief Returns the value of this SOGFunc_Impl evaluated at distance x.
	Real func( Real const x ) const;
	Real prob_sum_of_gaussians( Real const x ) const;

	/// @brief Returns the value of the first derivative of this SOGFunc_Impl at distance x.
	Real dfunc( Real const x ) const;

	void check_bounds( Real const x, Real const val ) const;
	void renormalize_weights();

	/// @brief show the definitio of this SOGFunc_Impl to the specified output stream.
	void show_definition( std::ostream & out ) const;

	/// @brief Initializes this SOGFunc_Impl from the given istream.
	/// @detailed The parameters are:
	/*!
	 * Initializes this SOGFunc_Impl from the given istream. An example
	 * of the type of string from which the istream should be constructed is:
	 * "SOGFUNC 2 19.396 7.643 0.4 17.312 2.4 0.6". The interpretation is to
	 * create initialize this SOGFunc_Impl object with the following parameters:
	 * - one Gaussian function with mean 19.396 and sd 7.643, weighted with a
	 * weight of 0.4
	 * - another Gaussian function with mean of 17.312 and sd 2.4, with a weight
	 * of 0.6.
	 * Weights need not add up to 1, but many times they will.
	 */
	void read_data( std::istream & in );

	void clear_(); // remove private data associated with this SOGFunc_Impl
	core::Real get_alt_score_( Real const x ) const;

	void upper_bound( const Real r );
	Real upper_bound() const;
	Real upper_bound_score() const;

	void smooth_to_zero( bool const setting );
	bool smooth_to_zero() const;

	void set_defaults();

	core::Real sog_cst_param() const;
	void sog_cst_param( core::Real const param );

private:
	utility::vector1< core::Real > means_;
	utility::vector1< core::Real > sdevs_;
	utility::vector1< core::Real > weights_;

	core::Real upper_bound_;
	core::Real score_upper_;

	core::Real sog_cst_param_;

	bool smooth_;
};

} // constraints
} // scoring
} // core

#endif
