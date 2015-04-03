// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/AtomTreeMultifunc.hh
/// @brief  Atom tree multifunction class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_AtomTreeMultifunc_hh
#define INCLUDED_core_optimization_AtomTreeMultifunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/optimization/MinimizerMap.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {

/// @brief Atom tree multifunction class
class AtomTreeMultifunc : public Multifunc {

public: // Creation

	// c-tor
	AtomTreeMultifunc(
		pose::Pose & pose_in,
		MinimizerMap & min_map_in,
		scoring::ScoreFunction const & scorefxn_in,
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	);

	/// @brief Destructor
	virtual ~AtomTreeMultifunc();

public: // Methods

	// func
	virtual
	Real
	operator ()( Multivec const & vars ) const;

	// dfunc
	virtual
	void
	dfunc( Multivec const & vars, Multivec & dE_dvars ) const;

	void set_deriv_check_result( NumericalDerivCheckResultOP deriv_check_result );

	/// @brief Error state reached -- derivative does not match gradient
	virtual
	void
	dump( Multivec const & vars, Multivec const & vars2 ) const;

protected: // accessors for subclasses
	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose() const;

	MinimizerMap const & min_map() const;

	scoring::ScoreFunction const & score_function() const;

private: // data

	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose_;

	/// non-const since min_map_ is modified by calls to dfunc()
	MinimizerMap & min_map_;

	scoring::ScoreFunction const & score_function_;

	bool deriv_check_;
	bool deriv_check_verbose_;
	NumericalDerivCheckResultOP deriv_check_result_;

}; // AtomTreeMultifunc


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_AtomTreeMultifunc_HH
