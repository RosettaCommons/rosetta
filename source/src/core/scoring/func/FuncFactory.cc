// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/func/FuncFactory.cc
/// @brief Factory for creating various types of constraints.
/// @author Greg Taylor <gktaylor@u.washington.edu>

// Unit headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>

// Package headers
#include <core/scoring/func/SumFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/USOGFunc.hh>
#include <core/scoring/func/SoedingFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/SigmoidFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/MixtureFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/CountViolFunc.hh>
#include <core/scoring/func/SkipViolFunc.hh>
#include <core/scoring/func/SquareWellFunc.hh>
#include <core/scoring/func/SquareWell2Func.hh>
#include <core/scoring/func/GaussianFunc.hh>
#include <core/scoring/func/ConstantFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/SplineFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/func/KarplusFunc.hh>
#include <core/scoring/func/IdentityFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/AmberPeriodicFunc.hh>

#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

using namespace core::scoring::constraints;

namespace core {
namespace scoring {
namespace func {

FuncFactory::~FuncFactory() = default;

void FuncFactory::add_type( std::string type_name, FuncOP new_func ) {
	func_types_[ type_name ] = new_func;
}

FuncOP FuncFactory::new_func( std::string const& type ) const {
	auto iter = func_types_.find( type );
	if ( iter != func_types_.end() ) {
		return iter->second->clone();
	} else {
		utility_exit_with_message("FuncFactory: unknown constraint function type: " + type );
		return nullptr;
	}
}

// initialization of functions which this factory knows how to instantiate
FuncFactory::FuncFactory(void) {
	FuncFactory::add_type( "HARMONIC", utility::pointer::make_shared< HarmonicFunc >(0,0) );
	FuncFactory::add_type( "SIGMOID", utility::pointer::make_shared< SigmoidFunc >(0,1) );
	FuncFactory::add_type( "AMBERPERIODIC", utility::pointer::make_shared< AmberPeriodicFunc >(0.0, 0.5, 1.0) );
	FuncFactory::add_type( "CIRCULARHARMONIC", utility::pointer::make_shared< CircularHarmonicFunc >(0,0) );
	FuncFactory::add_type( "MIXTUREFUNC", utility::pointer::make_shared< MixtureFunc >(0,0,0,0,0,0) );
	FuncFactory::add_type( "SCALARWEIGHTEDFUNC", utility::pointer::make_shared< ScalarWeightedFunc >(0,nullptr) );
	FuncFactory::add_type( "COUNTVIOLFUNC", utility::pointer::make_shared< CountViolFunc >(0,nullptr) );
	FuncFactory::add_type( "SKIPVIOLFUNC", utility::pointer::make_shared< SkipViolFunc >(0,nullptr) );
	FuncFactory::add_type( "GAUSSIANFUNC", utility::pointer::make_shared< GaussianFunc >(0,0) );
	FuncFactory::add_type( "CONSTANTFUNC", utility::pointer::make_shared< ConstantFunc >(0.) );
	FuncFactory::add_type( "BOUNDED", utility::pointer::make_shared< BoundFunc >(0,0,0,"dummy") );
	FuncFactory::add_type( "PERIODICBOUNDED", utility::pointer::make_shared< PeriodicBoundFunc >(0,0,0,"dummy",6.28) );
	FuncFactory::add_type( "OFFSETPERIODICBOUNDED", utility::pointer::make_shared< OffsetPeriodicBoundFunc >(0,0,0,"dummy",6.28,0.0) );
	FuncFactory::add_type( "SUMFUNC", utility::pointer::make_shared< SumFunc >() );
	FuncFactory::add_type( "SOGFUNC", utility::pointer::make_shared< SOGFunc >() );
	FuncFactory::add_type( "USOGFUNC", utility::pointer::make_shared< USOGFunc >() );
	FuncFactory::add_type( "SOEDINGFUNC", utility::pointer::make_shared< SoedingFunc >() );
	FuncFactory::add_type( "SPLINE", utility::pointer::make_shared< SplineFunc >() );
	FuncFactory::add_type( "SQUARE_WELL", utility::pointer::make_shared< SquareWellFunc >(0,0) );
	FuncFactory::add_type( "SQUARE_WELL2", utility::pointer::make_shared< SquareWell2Func >(0,0,0) );
	FuncFactory::add_type( "FADE", utility::pointer::make_shared< FadeFunc >(0,0,0) );
	FuncFactory::add_type( "LINEAR_PENALTY", utility::pointer::make_shared< LinearPenaltyFunction >(0,0,0,0) );
	FuncFactory::add_type( "KARPLUS", utility::pointer::make_shared< KarplusFunc >(6.98,-1.38,1.72,-1.05,0,0,0));
	FuncFactory::add_type( "IDENTITY", utility::pointer::make_shared< IdentityFunc >() );
	FuncFactory::add_type( "FLAT_HARMONIC", utility::pointer::make_shared< FlatHarmonicFunc >( 0, 0, 0 ) );
	FuncFactory::add_type( "TOPOUT", utility::pointer::make_shared< TopOutFunc >( 0, 0, 0 ) );
}

} //constraints
} //scoring
} //core
