// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/FitSimpleHelixMultiFunc.hh
/// @brief  Multifunction class implementation for fitting a simple (straight) helix to the Crick parameters.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_FitSimpleHelixMultiFunc_hh
#define INCLUDED_protocols_helical_bundle_FitSimpleHelixMultiFunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

//C++ headers
#include <string>

namespace protocols {
namespace helical_bundle {

/// @brief Multifunction class for fitting a simple (straight) helix to the Crick parameters.
///
class FitSimpleHelixMultiFunc : public core::optimization::Multifunc
{
public: // Creation

	/// @brief Destructor
	///
	~FitSimpleHelixMultiFunc() override;

	/// @brief Constructor
	///
	FitSimpleHelixMultiFunc(
		core::pose::Pose const &pose,
		std::string const & atom_name,
		core::Size const first_res_index,
		core::Size const res_per_repeat,
		core::Size const start_index,
		core::Size const end_index,
		core::Size const minimization_mode,
		core::Real const &rms_offset
	);

public: // Methods

	/// @brief Calculate function value (rms squared)
	///

	core::Real
	operator ()( core::optimization::Multivec const & vars ) const override;

	/// @brief Calculate function value derivatives
	///

	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const override;

	/// @brief Error state reached -- derivative does not match gradient
	///

	void
	dump( core::optimization::Multivec const & vars, core::optimization::Multivec const & vars2 ) const override;

private:

	/*************************************************************
	PRIVATE DATA
	*************************************************************/

	/// @brief A copy of the pose on which we'll be operating.  X,Y,Z coordinates of relevant residues will be extracted from this.
	///
	core::pose::Pose pose_;

	/// @brief The atoms that trace out the helix (e.g. "CA" or "N" or whatever)
	///
	std::string atom_name_;

	/// @brief The Rosetta-numbering index of the first residue in the helix that will be fit.
	///
	core::Size first_res_index_;

	/// @brief The number of residues per repeating unit in the helix.
	///
	core::Size residues_per_repeat_;

	/// @brief The first residue of the helix
	///
	core::Size start_index_;

	/// @brief The last residue of the helix
	///
	core::Size end_index_;

	/// @brief The minimization mode.
	/// @details If set to 0, the minimizer fits r1, omega1, and dz1, aligning the generated
	/// helix to the input pose at every step.  If set to 1, the minimizer fits r1, delta-omega1,
	/// and delta-z, with no alignment (i.e. the input pose must already be aligned to an ideal
	/// helix).
	core::Size minimization_mode_;

	/// @brief The small offset used by the RMS code to avoid zero determinants.
	/// @details This may need to be increased slightly to avoid bad geometry.
	core::Real rms_offset_;


	/*************************************************************
	PRIVATE MEMBER FUNCTIONS
	*************************************************************/

	/// @brief Convert the vars Multivec to the Crick parameters.
	///
	void vars_to_params( core::optimization::Multivec const &vars, core::Real &r1, core::Real &omega1, core::Real &dz1, core::Real &delta_omega1, core::Real &delta_z1 ) const;

	/// @brief Convert the Crick parameter derivatives to the derivative Multivec.
	///
	void params_derivs_to_vars( core::optimization::Multivec &deriv_vars, core::Real const &dE_dr1, core::Real const &dE_domega1, core::Real const &dE_ddz1, core::Real const &dE_ddelta_omega1, core::Real const &dE_ddelta_z1 ) const;


}; // FitSimpleHelixMultiFunc


} // namespace helical_bundle
} // namespace protocols


#endif // INCLUDED_protocols_helical_bundle_FitSimpleHelixMultiFunc_HH
