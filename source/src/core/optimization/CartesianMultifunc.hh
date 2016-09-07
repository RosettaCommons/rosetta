// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/CartesianMultifunc.hh
/// @brief
/// @author Frank DiMaio


#ifndef INCLUDED_core_optimization_CartesianMultifunc_hh
#define INCLUDED_core_optimization_CartesianMultifunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/optimization/CartesianMinimizerMap.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {

/// @brief Atom tree multifunction class
class CartesianMultifunc : public Multifunc {

public: // Creation

	// c-tor
	CartesianMultifunc(
		pose::Pose & pose_in,
		CartesianMinimizerMap & min_map_in,
		scoring::ScoreFunction const & scorefxn_in,
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	);

	/// @brief Destructor
	~CartesianMultifunc() override;

public: // Methods

	// func
	
	Real
	operator ()( Multivec const & vars ) const override;

	// dfunc
	
	void
	dfunc( Multivec const & vars, Multivec & dE_dvars ) const override;

	void set_deriv_check_result( NumericalDerivCheckResultOP deriv_check_result );

	/// @brief Error state reached -- derivative does not match gradient
	
	void
	dump( Multivec const & vars, Multivec const & vars2 ) const override;

protected: // accessors for subclasses
	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose() const;

	CartesianMinimizerMap const & min_map() const;

	scoring::ScoreFunction const & score_function() const;

private: // data

	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose_;

	/// non-const since min_map_ is modified by calls to dfunc()
	CartesianMinimizerMap & min_map_;

	scoring::ScoreFunction const & score_function_;

	bool deriv_check_;
	bool deriv_check_verbose_;
	NumericalDerivCheckResultOP deriv_check_result_;

}; // CartesianMultifunc


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_CartesianMultifunc_HH
