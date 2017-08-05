// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SmoothStepFunc.hh
/// @brief Definition for the smooth step function in the hydrate/SPaDES protocol.
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_core_scoring_func_SmoothStepFunc_hh
#define INCLUDED_core_scoring_func_SmoothStepFunc_hh

//#include <core/scoring/func/SmoothStepFunc.fwd.hh>

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


// C++ Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class SmoothStepFunc : public Func {
public:
	// The smooth step function is defined by a low and a high cut parameters, it will return:
	// if x < low; 0
	// if x > high; 1
	// if low < x < high; smooth step transition between 0 and 1
	SmoothStepFunc( Real const low, Real const high ): low_( low ), high_( high ){}

	/// @brief returns a clone of this SmoothStepFunc
	FuncOP clone() const { return FuncOP( new SmoothStepFunc( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

private:
	Real low_;
	Real high_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SmoothStepFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};



} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_SmoothStepFunc )
#endif // SERIALIZATION


#endif
