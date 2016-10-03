// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/func/CircularSplineFunc.hh
/// @brief  Similar to spline func but periodic on [0,360)
/// @author fpd

#ifndef INCLUDED_core_scoring_func_CircularSplineFunc_hh
#define INCLUDED_core_scoring_func_CircularSplineFunc_hh

// Unit Headers
#include <core/scoring/func/CircularSplineFunc.fwd.hh>

// Package Headers
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/Func.hh>

// Project Headers

// Utility and Numeric Headers
#include <utility/vector1.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>

// C++ Headers
#include <iostream>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class CircularSplineFunc : public Func {
public:

	CircularSplineFunc() : weight_(0.0), convert_to_degrees_( false ) {}
	CircularSplineFunc( core::Real weight_in, utility::vector1< core::Real> energies_in,
		bool convert_to_degrees = false ); // true might be better? esp. if used in DihedralConstraint.

	~CircularSplineFunc() {}

	virtual
	FuncOP clone() const { return FuncOP( new CircularSplineFunc( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	core::Real get_weight() { return weight_; }

	void train( utility::vector1< core::Real> energies_in );

	virtual
	void read_data ( std::istream &in );

	virtual
	core::Real func( core::Real const x ) const;

	virtual
	core::Real dfunc( core::Real const x ) const;

	virtual
	void show_definition( std::ostream &out ) const;

	virtual
	core::Size show_violations( std::ostream &out, core::Real x, core::Size verbose_level, core::Real threshold = 1 ) const;

private:

	core::Real weight_;
	bool convert_to_degrees_;
	numeric::interpolation::spline::CubicSpline spline_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // CircularSplineFunc class

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_CircularSplineFunc )
#endif // SERIALIZATION


#endif
