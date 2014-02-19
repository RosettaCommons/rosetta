// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/normalmode/NormalModeMultifunc.hh
/// @brief  NormalMode multifunction class
/// @author Hahnbeom Park


#ifndef INCLUDED_protocols_normalmode_NormalModeMultifunc_hh
#define INCLUDED_protocols_normalmode_NormalModeMultifunc_hh

// Package headers
#include <protocols/normalmode/NormalMode.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/optimization/MinimizerMap.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace normalmode {

using namespace core;
using namespace core::optimization;

/// @brief Atom tree multifunction class
class NormalModeMultifunc : public Multifunc
{
public: // Creation

	/// @brief Destructor
	virtual ~NormalModeMultifunc();

	// c-tor
	NormalModeMultifunc(
		pose::Pose & pose_in,
		MinimizerMap & min_map_in,
		scoring::ScoreFunction const & scorefxn_in,
		protocols::normalmode::NormalMode const & normalmode_in,
		bool const use_omega = true,
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	);

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

	void set_modes( utility::vector1< Size > );

	Multivec 
	dofs_to_vars( Multivec const & dofs ) const;

	Multivec 
	vars_to_dofs( Multivec const & vars ) const;

	Size nmodes() const { return modes_using_.size(); }

	Size nvar() const { return nvar_; }

private:

	void get_dofs_for_pose0( );
	void get_dofs_map( );

	void set_default_modes();

	Real dofs_for_pose0( Size const i_dof ) const { return dofs_for_pose0_[i_dof]; }

	Multivec
	dEddofs_to_dEdvars( Multivec const & dEdtors ) const;

	Real get_modescale( Size const modeno ) const;

protected: // accessors for subclasses
	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose() const;

	MinimizerMap const & min_map() const;

	scoring::ScoreFunction const & score_function() const;

private: // data

	// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose_;

	// non-const since min_map_ is modified by calls to dfunc()
	MinimizerMap & min_map_;

	scoring::ScoreFunction const & score_function_;

	// Number of variables
	Size nvar_;

	// Whether to use omega as variable?
	bool use_omega_;

	// Dampening factor for normal modes
	Real k_dampen_;

	// Reference pose during minimization( set as initial structure )
	pose::Pose & pose0_;
	Multivec dofs_for_pose0_;

	// Normalmode
	protocols::normalmode::NormalMode const NM_;
	utility::vector1< Size > modes_using_;

	// Map between NM torsionID <-> min_map DOF_ID
	std::map< Size, Size > map_NM_to_DOF_;
	std::map< Size, Size > map_DOF_to_NM_;

	// Map between vars ID <-> min_map DOF_ID
	std::map< Size, Size > map_var_to_DOF_;

	// Map between vars ID <-> NormalMode number
	std::map< Size, Size > map_var_to_modeno_;

	// vars type
	utility::vector1< std::string > var_type_;

	bool deriv_check_;
	bool deriv_check_verbose_;
	NumericalDerivCheckResultOP deriv_check_result_;

}; // NormalModeMultifunc


} // namespace normalmode
} // namespace protocols


#endif // INCLUDED_protocols_normalmode_NormalModeMultifunc_HH
