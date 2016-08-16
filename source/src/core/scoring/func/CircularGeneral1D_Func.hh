// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/CircularGeneral1D_Func.hh
/// @brief A general 1D function that can be initialized by FArray or from text file. There's stuff like this all over Rosetta.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_func_CircularGeneral1D_Func_HH
#define INCLUDED_core_scoring_func_CircularGeneral1D_Func_HH

#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Function that allows return of arbitrary FArrays -- this time circularized.
class CircularGeneral1D_Func : public Func {
public:
	CircularGeneral1D_Func(
		ObjexxFCL::FArray1D< core::Real > const & data,
		core::Real const & xmin,
		core::Real const & xbin
	);

	CircularGeneral1D_Func( std::string const & filename );

	FuncOP clone() const { return FuncOP( new CircularGeneral1D_Func( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

private:

	ObjexxFCL::FArray1D< core::Real > data_;
	Real xmin_, xbin_;
	Size num_bins_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CircularGeneral1D_Func();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_CircularGeneral1D_Func )
#endif // SERIALIZATION


#endif
