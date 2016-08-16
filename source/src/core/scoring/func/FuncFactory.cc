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

FuncFactory::~FuncFactory() {}

void FuncFactory::add_type( std::string type_name, FuncOP new_func ) {
	func_types_[ type_name ] = new_func;
}

FuncOP FuncFactory::new_func( std::string const& type ) const {
	FuncTypes::const_iterator iter = func_types_.find( type );
	if ( iter != func_types_.end() ) {
		return iter->second->clone();
	} else {
		utility_exit_with_message("FuncFactory: unknown constraint function type: " + type );
		return NULL;
	}
}

// initialization of functions which this factory knows how to instantiate
FuncFactory::FuncFactory(void) {
	FuncFactory::add_type( "HARMONIC", FuncOP( new HarmonicFunc(0,0) ) );
	FuncFactory::add_type( "SIGMOID", FuncOP( new SigmoidFunc(0,1) ) );
	FuncFactory::add_type( "AMBERPERIODIC", FuncOP( new AmberPeriodicFunc(0.0, 0.5, 1.0) ) );
	FuncFactory::add_type( "CIRCULARHARMONIC", FuncOP( new CircularHarmonicFunc(0,0) ) );
	FuncFactory::add_type( "MIXTUREFUNC", FuncOP( new MixtureFunc(0,0,0,0,0,0) ) );
	FuncFactory::add_type( "SCALARWEIGHTEDFUNC", FuncOP( new ScalarWeightedFunc(0,0) ) );
	FuncFactory::add_type( "COUNTVIOLFUNC", FuncOP( new CountViolFunc(0,0) ) );
	FuncFactory::add_type( "SKIPVIOLFUNC", FuncOP( new SkipViolFunc(0,0) ) );
	FuncFactory::add_type( "GAUSSIANFUNC", FuncOP( new GaussianFunc(0,0) ) );
	FuncFactory::add_type( "CONSTANTFUNC", FuncOP( new ConstantFunc(0.) ) );
	FuncFactory::add_type( "BOUNDED", FuncOP( new BoundFunc(0,0,0,"dummy") ) );
	FuncFactory::add_type( "PERIODICBOUNDED", FuncOP( new PeriodicBoundFunc(0,0,0,"dummy",6.28) ) );
	FuncFactory::add_type( "OFFSETPERIODICBOUNDED", FuncOP( new OffsetPeriodicBoundFunc(0,0,0,"dummy",6.28,0.0) ) );
	FuncFactory::add_type( "SUMFUNC", FuncOP( new SumFunc() ) );
	FuncFactory::add_type( "SOGFUNC", FuncOP( new SOGFunc() ) );
	FuncFactory::add_type( "USOGFUNC", FuncOP( new USOGFunc() ) );
	FuncFactory::add_type( "SOEDINGFUNC", FuncOP( new SoedingFunc() ) );
	FuncFactory::add_type( "SPLINE", FuncOP( new SplineFunc() ) );
	FuncFactory::add_type( "SQUARE_WELL", FuncOP( new SquareWellFunc(0,0) ) );
	FuncFactory::add_type( "SQUARE_WELL2", FuncOP( new SquareWell2Func(0,0,0) ) );
	FuncFactory::add_type( "FADE", FuncOP( new FadeFunc(0,0,0) ) );
	FuncFactory::add_type( "LINEAR_PENALTY", FuncOP( new LinearPenaltyFunction(0,0,0,0) ) );
	FuncFactory::add_type( "KARPLUS", FuncOP( new KarplusFunc(6.98,-1.38,1.72,-1.05,0,0,0) ));
	FuncFactory::add_type( "IDENTITY", FuncOP( new IdentityFunc() ) );
	FuncFactory::add_type( "FLAT_HARMONIC", FuncOP( new FlatHarmonicFunc( 0, 0, 0 ) ) );
	FuncFactory::add_type( "TOPOUT", FuncOP( new TopOutFunc( 0, 0, 0 ) ) );
}

} //constraints
} //scoring
} //core
