// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/TopOutFunc.hh
/// @brief Implementation of phenix "top-out" function
///   Similar to Geman-McClure: harmonic near 'x0_', flat past 'limit_'
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_func_TopOutFunc_hh
#define INCLUDED_core_scoring_func_TopOutFunc_hh

#include <core/scoring/func/TopOutFunc.fwd.hh>
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

class TopOutFunc : public Func {
public:
	TopOutFunc( Real weight_in, Real x0_in, Real limit_in );

	virtual FuncOP clone() const;
	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	virtual Real func( Real const x ) const;
	virtual Real dfunc( Real const x ) const;

	virtual void read_data( std::istream & in );

	virtual void show_definition( std::ostream &out ) const;

	Real x0() const { return x0_; }
	Real limit() const { return limit_; }
	Real weight() const { return weight_; }

	void x0( Real x ) { x0_ = x; }
	void limit( Real limit ) { limit_ = limit; }
	void weight( Real weight ) { weight_ = weight; }

	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	Real x0_;
	Real weight_;
	Real limit_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	TopOutFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_TopOutFunc )
#endif // SERIALIZATION


#endif
