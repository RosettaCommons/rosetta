// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularPowerFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_func_CircularPowerFunc_hh
#define INCLUDED_core_scoring_func_CircularPowerFunc_hh

#include <core/scoring/func/Func.hh>
#include <core/types.hh>

// C++ Headers
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Generalization of CircularCircularPowerFunc -- other exponents allowed.
/// @brief Operates in radians, like DihedralConstraint.
class CircularPowerFunc : public Func {
public:
	CircularPowerFunc( Real const x0_radians, Real const sd_radians, int const power, Real const weight );

	FuncOP
	clone() const;

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real
	func( Real const x ) const;

	Real
	dfunc( Real const x ) const;

private:
	Real const x0_;
	Real const sd_;
	int const power_;
	Real const weight_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< CircularPowerFunc > & construct );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_CircularPowerFunc )
#endif // SERIALIZATION


#endif
