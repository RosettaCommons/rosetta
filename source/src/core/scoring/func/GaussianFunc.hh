// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_func_GaussianFunc_hh
#define INCLUDED_core_scoring_func_GaussianFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Derived class of class Func representing a Gaussian distribution with a user-specified
/// mean and standard deviation.
class GaussianFunc : public Func {
public:

	/*!
	* Constuctor for GaussianFunc. Arguments to the constructor are:
	* - mean: parameter representing the mean of this function.
	* - sd: parameter representing the standard deviation of this function.
	*/

	GaussianFunc (
		Real const mean,
		Real const sd
	) :
		mean_         ( mean ),
		sd_           ( sd ),
		use_log_score_( true ),
		weight_ (1.0)
	{}

	/// @brief returns a clone of this GaussianFunc
	FuncOP clone() const { return FuncOP( new GaussianFunc( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	/// @brief Returns the value of this GaussianFunc evaluated at distance x.
	Real func( Real const x ) const;

	/// @brief Returns the value of the first derivative of this GaussianFunc at distance x.
	Real dfunc( Real const x ) const;

	/// @brief show the definitio of this GaussianFunc to the specified output stream.
	virtual void show_definition( std::ostream &out ) const;

	/// @brief Calls show( out ) on this GaussianFunc.
	friend std::ostream& operator<<(std::ostream& out, const GaussianFunc& f ) {
		f.show( out );
		return out;
	} // operator<<


	/// @brief
	/// The parameters are:
	/*!
	* Initializes this GaussianFunc from the given istream. An example
	* of the type of string from which the istream should be constructed is:
	* "GAUSSIANFUNC 19.396 7.643". The interpretation is to
	* create initialize this GaussianFunc object with the following parameters:
	* - mean 19.396
	* - sd 7.643
	*/
	void read_data( std::istream& in );

private:
	Real mean_;
	Real sd_;
	bool use_log_score_;

	/// @brief A multiplier for the Gaussian function.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	Real weight_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	GaussianFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_GaussianFunc )
#endif // SERIALIZATION


#endif
