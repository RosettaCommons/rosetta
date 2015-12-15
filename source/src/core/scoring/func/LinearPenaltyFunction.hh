// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/LinearPenaltyFunction.hh
/// @brief Linear penalty for values outside a certain range
/// @author Dominik Gront

#ifndef INCLUDED_core_scoring_func_LinearPenaltyFunction_hh
#define INCLUDED_core_scoring_func_LinearPenaltyFunction_hh

#include <core/scoring/func/Func.hh>

#include <core/types.hh>

// #include <utility/pointer/ReferenceCount.hh>

// #include <ObjexxFCL/format.hh>
// #include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

class LinearPenaltyFunction : public Func {
public:
	LinearPenaltyFunction( Real const x_middle, Real const well_depth, Real const half_width, Real const slope ):
		x_middle_( x_middle ), well_depth_( well_depth ), half_width_ ( half_width ), slope_ ( slope ) {}

	FuncOP
	clone() const { return FuncOP( new LinearPenaltyFunction( *this ) ); }

	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream& in );

	void show_definition( std::ostream &out ) const;

	Real get_x_middle() const { return x_middle_; }
	Real get_well_depth() const { return well_depth_; }
	Real get_half_width() const { return half_width_; }
	Real get_slope() const { return slope_; }
	void set_x_middle( Real x ) { x_middle_ = x; }
	void set_well_depth( Real well_depth ){ well_depth_ = well_depth;}
	void set_half_width(Real half_width) { half_width_ = half_width;}
	void set_slope(Real slope) { slope_ = slope; }
	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	Real x_middle_;
	Real well_depth_;
	Real half_width_;
	Real slope_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	LinearPenaltyFunction();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_LinearPenaltyFunction )
#endif // SERIALIZATION


#endif
