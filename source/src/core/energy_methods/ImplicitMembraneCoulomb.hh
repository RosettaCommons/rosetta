// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/ImplicitMembraneCoulomb.hh
/// @brief Minimal class for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)


#ifndef INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_hh
#define INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_hh

#include <core/energy_methods/ImplicitMembraneCoulomb.fwd.hh>
#include <core/types.hh>
#include <numeric/cubic_polynomial.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace core {
namespace energy_methods {

/// @brief Minimal class for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)
class ImplicitMembraneCoulomb : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	ImplicitMembraneCoulomb();

	/// @brief Copy constructor.
	ImplicitMembraneCoulomb(ImplicitMembraneCoulomb const & src);

	/// @brief Destructor.
	~ImplicitMembraneCoulomb() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	ImplicitMembraneCoulombOP clone() const;

	/// @brief Initialize max distance and coulomb constants
	void initialize();

	/// @brief Compute depth- and bilayer-dependent dielectric constant
	core::Real
	compute_depth_and_bilayer_dep_dielectric(
		core::Real const fi, /* fractional hydration of atom i */
		core::Real const fj, /* fractional hydration of atom j */
		core::Real const d /* Distance between atoms i and j */
	) const;

	///@brief Compute the distance dependent sigmoid function
	core::Real
	sigmoidal_eps( core::Real const sigmoid_D0,/*upper limit of sigmoid*/
		core::Real const sigmoid_D, /*lower limit of sigmoid*/
		core::Real const sigmoid_s,/*slope of sigmoid*/
		core::Real const d) const;

	core::Real
	sigmoidal_deps_dr( core::Real const sigmoid_D0,/*upper limit of sigmoid*/
		core::Real const sigmoid_D, /*lower limit of sigmoid*/
		core::Real const sigmoid_s,/*slope of sigmoid*/
		core::Real const d) const;

	/// @brief Compute the derivative of depth- and bilayer-dependent dielectric
	core::Real
	compute_deps_dr(
		core::Real const fi, /* fractional hydration of atom i */
		core::Real const fj, /* fractional hydration of atom j */
		core::Real const d
	) const;

	/// @brief Evaluate atom pair electrostatics energy (pass through)
	core::Real
	eval_atom_atom_fa_elecE(
		core::Vector const & i_xyz,
		core::Real const i_charge,
		core::Real const i_hyd,
		core::Vector const & j_xyz,
		core::Real const j_charge,
		core::Real const j_hyd
	) const;

	/// @brief Evaluate atom pair electrostatics energy
	core::Real
	eval_atom_atom_fa_elecE(
		core::Vector const & i_xyz,
		core::Real const i_charge,
		core::Real const i_hyd,
		core::Vector const & j_xyz,
		core::Real const j_charge,
		core::Real const j_hyd,
		DistanceSquared & d2
	) const;

	/// @brief Get the key numeric value for derivative calculations
	core::Real
	eval_dfa_elecE_dr_over_r(
		core::Real const dis2,
		core::Real const q1,
		core::Real const q2,
		core::Real const i_hyd,
		core::Real const j_hyd
	) const;

	core::Vector
	eval_dfa_elecE_df(
		core::Real const dis2,
		core::Real const q1,
		core::Real const q2,
		core::Real const i_hyd,
		core::Real const j_hyd,
		core::Vector const di_hyd_df,
		core::Vector const dj_hyd_df
	) const;

	numeric::CubicPolynomial compute_hipoly(
		core::Real const fi,
		core::Real const fj
	)const;

	numeric::CubicPolynomial compute_lowpoly(
		core::Real const fi,
		core::Real const fj
	)const;

	core::Real compute_min_dis_score(
		core::Real const fi,
		core::Real const fj
	)const;

private:

	core::Real max_dis_;
	core::Real max_dis2_;
	core::Real min_dis_;
	core::Real min_dis2_;

	core::Real low_poly_start_;
	core::Real low_poly_start2_;
	core::Real low_poly_end_;
	core::Real low_poly_end2_;

	core::Real hi_poly_start_;
	core::Real hi_poly_start2_;
	core::Real hi_poly_end_;
	core::Real hi_poly_end2_;

	core::Real C0_;
	core::Real C1_;
	core::Real dEfac_;

};

} //energy_methods
} //core

#endif //INCLUDED_core_energy_methods_ImplicitMembraneCoulomb_hh
