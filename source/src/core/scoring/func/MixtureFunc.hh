// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/MixtureFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_func_MixtureFunc_hh
#define INCLUDED_core_scoring_func_MixtureFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Derived class of class Func representing a Mixture of several distinct functions. The
/// function is of the form ln( g(r) / h(r) ), where g(r) is a mixture of a Gaussian and Exponential
/// distributions, and h(r) is a Gaussian distribution. See methods and implementation for more
/// information.
class MixtureFunc : public Func {
public:

	/*!
	* Constructor for MixtureFunc. Arguments to the constructor are:
	* - anchor: parameter representing the value at which this function is anchored, represents the
	* mean of the Gaussian distribution and the highest point of the exponential distribution.
	* - gaussian_param: parameter for Gaussian portion of g(r), representing the standard deviation
	* of a Gaussian distribution around anchor.
	* - exp_param: parameter for Exponential portion of g(r), representing the rate at which
	* the exponential distribution drops off from anchor.
	* - mixture_param: parameter describing the mixture of the Gaussian and Exponential functions
	*   that make up g(r).
	* - bg_mean: parameter representing the mean of h(r).
	* - bg_sd: parameter representing the standard deviation of h(r).
	*/

	MixtureFunc (
		Real const anchor,
		Real const gaussian_param,
		Real const exp_param,
		Real const mixture_param,
		Real const bg_mean,
		Real const bg_sd
	) :
		fmax_          ( 0.0 ), // this is always read, but not always initialized.  Is this a good value?
		anchor_        ( anchor ),
		gaussian_param_( gaussian_param ),
		exp_param_     ( exp_param ),
		mixture_param_ ( mixture_param ),
		bg_mean_       ( bg_mean ),
		bg_sd_         ( bg_sd )
	{
		if ( anchor_ > 1e-10 ) {
			verify_parameters_();
		}
	}

	/// @brief returns a clone of this MixtureFunc
	FuncOP clone() const { return FuncOP( new MixtureFunc( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	/// @brief Returns the value of this MixtureFunc evaluated at distance x.
	Real func( Real const x ) const;

	/// @brief Returns the value of the first derivative of this MixtureFunc at distance x.
	Real dfunc( Real const x ) const;

	/// @brief show the definitio of this MixtureFunc to the specified output stream.
	virtual void show_definition( std::ostream &out ) const;

	/// @brief Function that's used for debugging. Given x, this calculates
	/// g(x), h(x), g'(x) and h'(x).
	Real
	dfunc_component(
		Real const x,
		Real & g,
		Real & h,
		Real & g_prime,
		Real & h_prime
	) const;

	/// @brief Calculates the K-L divergence between the inferred and background distributions.
	Real calc_kl_divergence() const;

	/// @brief Prints this MixtureFunc to the given ostream.
	virtual void show( std::ostream& out ) const;

	/// @brief Calls show( out ) on this MixtureFunc.
	friend std::ostream& operator<<( std::ostream& out, const MixtureFunc& f ) {
		f.show( out );
		return out;
	} // operator<<

	/// @brief
	/// The parameters are:
	/*!
	* Initializes this MixtureFunc from the given istream. An example
	* of the type of string from which the istream should be constructed is:
	* "MIXTUREFUNC 6.9734 3.598 0.222 0.872 19.396 7.643". The interpretation is to
	* create initialize this MixtureFunc object with the following parameters:
	* - anchor 6.9734
	* - gaussian_param 3.598
	* - exp_param 0.222
	* - mixture_param 0.872
	* - bg_mean 19.396
	* - bg_sd 7.643
	*/
	void read_data( std::istream& in );

	/// @brief Returns the value of this MixtureFunc evaluated at distance x.
	Real func_( Real x ) const;


private:
	void verify_parameters_();

	// Real distance_cutoff_;
	Real rmax_;
	Real fmax_;
	Real anchor_;
	Real gaussian_param_;
	Real exp_param_;
	Real mixture_param_;
	Real bg_mean_;
	Real bg_sd_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	MixtureFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_MixtureFunc )
#endif // SERIALIZATION


#endif
