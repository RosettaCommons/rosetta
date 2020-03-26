// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/CharmmPeriodicFunc.hh
/// @brief Definition for a periodic function used in constraints.
/// @brief feel free to add definitions for other (derived) periodic functions to this file
/// @author Florian Richter

#ifndef INCLUDED_core_scoring_func_CharmmPeriodicFunc_hh
#define INCLUDED_core_scoring_func_CharmmPeriodicFunc_hh

#include <core/scoring/func/Func.hh>

#include <utility/VirtualBase.hh>
#include <utility/pointer/owning_ptr.hh>

//#include <numeric/angle.functions.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {


/// @brief function of type y = 0.5 * k * (1 - cos(n * (x - x0) ) )
class CharmmPeriodicFunc : public Func {

public:
	CharmmPeriodicFunc ( Real const x0_in, Real const k_in, Real const n_periodic_in) : x0_( x0_in ), k_( k_in ), n_periodic_( n_periodic_in ){}

	FuncOP
	clone() const override { return utility::pointer::make_shared< CharmmPeriodicFunc >( *this ); }

	bool operator == ( Func const & other ) const override;
	bool same_type_as_me( Func const & other ) const override;

	Real func( Real const x ) const override;
	Real dfunc( Real const x ) const override;

	void read_data( std::istream& in ) override;

	void show_definition( std::ostream &out ) const override;

	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const override;


private:
	Real x0_;
	Real k_;
	Real n_periodic_;


#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CharmmPeriodicFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //charmm periodic func


} // constraints
} // scoring
} // core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_CharmmPeriodicFunc )
#endif // SERIALIZATION


#endif
