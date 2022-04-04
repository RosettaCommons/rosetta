// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/lkball/LK_DomeEnergy.hh
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_scoring_lkball_LK_DomeEnergy_hh
#define INCLUDED_core_scoring_lkball_LK_DomeEnergy_hh

// Unit headers
#include <core/scoring/lkball/LK_DomeEnergy.fwd.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/lkball/LK_DomeInfo.hh>

// Package headers
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/ResidueType.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <numeric/cubic_polynomial.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


// forward declaration for testing
class LK_DomeEnergyTests;



namespace core {
namespace scoring {
namespace lkball {


class LK_DomeEnergy : public scoring::methods::ContextDependentTwoBodyEnergy {

	friend class ::LK_DomeEnergyTests;

public:
	LK_DomeEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	scoring::methods::EnergyMethodOP
	clone() const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	Distance          // Just assume the worst that max_angle is 180
	// lk_dome  Atom 2.65 LK_ball 2.65 LK_dome w_rad(1.4) fade(1.92) radius(2.2) Atom
	// lk_dome = 10.82
	// lk_dome_br  Atom 2.65 LK_ball 2.65 LK_dome fade(2.23) LK_ball 2.65 Atom
	// lk_dome_br = 10.18
	atomic_interaction_cutoff() const override;
	// water_water_fade = 2.23
	// lk_Dome_bridge = 5.3 + 2.23 + 2.65 = 10.18

	Real                                            // only true if the above assumes max_angle is 180
	water_atom_interaction_cutoff() const;

	Real                                            // only true if the above assumes max_angle is 180
	water_atom_interaction_cutoff2() const;

	Real                                            // only true if the above assumes max_angle is 180
	water_water_bridge_cutoff2() const;



	Real
	get_sol_value( chemical::AtomType const & at, Size nattached_waters ) const;

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const override;


	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const override;


	bool
	defines_intrares_energy( EnergyMap const & weights ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	core::Size version() const override;

	// We don't want any of these
	void indicate_required_context_graphs(
		utility::vector1< bool > &
	) const override {}

	Real
	dome_water_dist2(
		Vector const & base_atom_xyz,
		Vector const & water,
		Vector const & other_atom_xyz,
		Real & water_atom_distance,
		Vector & wadj_water
	) const;

	Vector
	get_dome_water(
		Vector const & base_atom_xyz,
		Vector const & water,
		Vector const & other_atom_xyz,
		Real & water_atom_distance,
		Vector & wadj_water
	) const;

	Real
	eval_lk_fraction( Real const d2_delta, Real const width ) const;

	Real
	get_lkd_frac_dome_dist( Real dw_dist, Real atom2_type ) const;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		LKD_ResidueInfo const & rsd1_info,
		conformation::Residue const & rsd2,
		LKD_ResidueInfo const & rsd2_info,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	void
	setup_d2_bounds();

	void
	accumulate_single_atom_contributions(
		Size const atom1,
		Size const,
		Size atom1_n_attached_waters,
		Size atom1_start_water,
		WaterCoords const & atom1_waters,
		WaterOcclusions const & atom1_occlusions,
		conformation::Residue const & rsd1,
		Size const atom2_type_index,
		Vector const & atom2_xyz,
		Real const lk_desolvation_of_atom1_by_water,
		EnergyMap & emap
	) const;


	bool
	requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose::Pose const & ) const override;

	void
	my_setup_for_scoring_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const &,
		ScoreFunction const &,
		basic::datacache::BasicDataCache & residue_data_cache
	) const;

	bool
	requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const override;

	void
	setup_for_derivatives_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const &,
		ScoreFunction const &,
		ResSingleMinimizationData &,
		basic::datacache::BasicDataCache & residue_data_cache
	) const override;

	void
	setup_for_scoring(
		pose::Pose &,
		scoring::ScoreFunction const &
	) const override;


	void
	setup_for_packing( pose::Pose &, utility::vector1< bool > const &, utility::vector1< bool > const & ) const override;

	void
	setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const &) const override;


	void
	finalize_after_minimizing(
		pose::Pose & pose
	) const override;


	void
	fill_static_occlusions( pose::Pose & pose ) const;

	void
	fill_occlusions_1way(
		conformation::Residue const & rsd1,
		LKD_ResidueInfo & rsd1_info,
		conformation::Residue const & rsd2
	) const;

	Real
	get_lkd_bridge_fractional_1way(
		Vector const & atom1_xyz,
		Size n_atom1_waters,
		Size atom1_start_water,
		Size n_atom2_waters,
		Size atom2_start_water,
		WaterCoords const & atom1_waters,
		WaterCoords const & atom2_waters,
		WaterOcclusions const & atom1_occlusions
	) const;

	void
	single_water_atom_fractions(
		Vector const & base_atom_xyz,
		Vector const & water_xyz,
		Vector const & other_atom_xyz,
		Size other_atom_type_index,
		Real & lk_dome_iso_frac,
		Real & lk_dome_frac
	) const;

	Real
	single_water_water_fraction_1way(
		Vector const & base_atom_xyz,
		Vector const & water_xyz,
		Vector const & other_water_xyz
	) const;


	Real
	get_lkb_bridge2_fractional_1way(
		Size n_atom1_waters,
		Size atom1_start_water,
		Size n_atom2_waters,
		Size atom2_start_water,
		WaterCoords const & atom1_waters,
		WaterCoords const & atom2_waters,
		WaterOcclusions const & atom1_occlusions
	) const;

	Real
	single_water_water_bridge2_fraction_1way(
		Vector const & water_xyz,
		Vector const & other_water_xyz
	) const;

	void
	update_cached_lkb_resinfo(
		conformation::Residue const & rsd,
		basic::datacache::BasicDataCache & residue_data_cache,
		bool compute_derivs,
		bool clear_occlusions
	) const;

	LK_BallEnergyCOP lk_ball() const;

	Real occlusion_max() const;
	Real occlusion_min() const;

	Real w_dist() const;

	Real water_adjust() const;


	Real get_avail( Real occl ) const;

	Real
	dlk_frac_sigmoid( Real distance ) const;

	Real
	dmy_lk_fraction( Real distance ) const;

	Real
	lk_frac_sigmoid( Real distance ) const;

	Real
	my_lk_fraction( Real distance ) const;

	void
	setup_poly_params();

	void
	sum_deriv_contributions_for_heavyatom_pair_one_way(
		Size const heavyatom1,
		conformation::Residue const & rsd1,
		LKD_ResidueInfo const & rsd1_info,
		Size const heavyatom2,
		conformation::Residue const & rsd2,
		LKD_ResidueInfo const & rsd2_info,
		EnergyMap const & weights,
		Real const weight_factor,
		Real const d2,
		utility::vector1< DerivVectorPair > & r1_at_derivs,
		utility::vector1< DerivVectorPair > & r2_at_derivs
	) const;

	void
	sum_deriv_contributions_for_heavyatom_pair(
		Real const d2,
		Size const heavyatom1,
		conformation::Residue const & rsd1,
		LKD_ResidueInfo const & rsd1_info,
		Size const heavyatom2,
		conformation::Residue const & rsd2,
		LKD_ResidueInfo const & rsd2_info,
		pose::Pose const &,
		EnergyMap const & weights,
		Real const cp_weight,
		utility::vector1< DerivVectorPair > & r1_at_derivs,
		utility::vector1< DerivVectorPair > & r2_at_derivs
	) const;

	Real
	eval_d_lk_fraction_dr_over_r( Real const d2_delta, Real const width ) const;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & , //pose,
		ScoreFunction const & , //scorefxn,
		kinematics::MinimizerMapBase const &, // min_map,
		ResSingleMinimizationData const &,// res1data,
		ResSingleMinimizationData const &,// res2data,
		ResPairMinimizationData & pair_data
	) const override;

	mutable bool dump_dome_waters_;
	mutable bool dump_dome_bridge_waters_;
	mutable bool debug_disable_count_pair_;

private:

	mutable bool packing_;
	mutable bool minimizing_;

	Real h2o_radius_;
	Real ramp_width_A2_;
	Real overlap_width_A2_;
	Real ball_overlap_width_A2_;
	Real w_dist_;
	Real water_adjust_;

	LK_BallEnergyCOP lk_ball_;
	Real max_angle_cos_;
	Real min_angle_cos_;
	Real max_angle_sin_;
	Real min_angle_sin_;
	Real occlusion_max_;
	Real occlusion_min_;

	numeric::CubicPolynomial poly_params_;


	utility::vector1< Real > d2_low_;

	mutable utility::vector1<numeric::xyzVector<Real>> debug_dome_waters_;

};


void
evaluate_lk_dome_energy_for_atom_ranges(
	LK_DomeEnergy const & lk_dome,
	conformation::Residue const & rsd1,
	LKD_ResidueInfo const & rsd1_info,
	conformation::Residue const & rsd2,
	LKD_ResidueInfo const & rsd2_info,
	ScoreFunction const & sfxn,
	etable::count_pair::CountPairFunction const & cpfxn,
	Size const res1_start_atom,
	Size const res1_end_atom,
	Size const res2_start_atom,
	Size const res2_end_atom,
	EnergyMap & emap
);


struct DerivativeFinder {

	Vector B;
	Vector W;
	Vector Ot;

	Vector W_to_Ot;
	Real W_to_Ot_norm;
	Real W_to_Ot_norm2;
	Real W_to_Ot_norm3;
	Vector B_to_W;
	Real B_to_W_norm;
	Real B_to_W_norm2;
	Real B_to_W_norm3;
	Real norm_prod;

	Real cose;
	Real sine;

	Real w_dist_;
	Real wadj;

	bool colinear_water;

	utility::vector1<Real> x;


	DerivativeFinder( Vector const & base, Vector const & water, Vector const & other,
		Real min_cos, Real min_sin, Real max_cos, Real max_sine, Real w_dist, Real water_adjust );

	numeric::xyzMatrix<Real>
	dDome_dBase();

	numeric::xyzMatrix<Real>
	dDome_dWater();

	numeric::xyzMatrix<Real>
	dDome_dOther();

	numeric::xyzMatrix<Real>
	colinear_dDome_dBase();

	numeric::xyzMatrix<Real>
	colinear_dDome_dWater();

	numeric::xyzMatrix<Real>
	colinear_dDome_dOther();

	numeric::xyzMatrix<Real>
	fringe_dDome_dBase();

	numeric::xyzMatrix<Real>
	fringe_dDome_dWater();

	numeric::xyzMatrix<Real>
	fringe_dDome_dOther();
};



struct DerivativeFinderWadj {

	Vector B;
	Vector W;

	Vector B_to_W;
	Real B_to_W_norm;
	Real B_to_W_norm2;
	Real B_to_W_norm3;
	Real wadj;
	Real w_dist_;


	utility::vector1<Real> x;


	DerivativeFinderWadj( Vector const & base, Vector const & water, Real w_dist, Real water_adjust );

	numeric::xyzMatrix<Real>
	dWadj_dBase();

	numeric::xyzMatrix<Real>
	dWadj_dWater();

};


} // guidance_scoreterms
} // pack
} // core

#endif
