// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/LK_BallEnergy.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_LK_BALLENERGY_HH
#define INCLUDED_core_scoring_methods_LK_BALLENERGY_HH

// Unit Headers
#include <core/scoring/methods/LK_BallEnergy.fwd.hh>
#include <core/scoring/methods/LK_BallInfo.hh>

#include <core/scoring/methods/GenBornEnergy.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/hbonds/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {



///
class LK_BallEnergy : public ContextIndependentTwoBodyEnergy {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:
	/// convenience typedefs
	typedef chemical::ResidueType ResidueType;
	typedef utility::vector1< Size > Sizes;
	typedef utility::vector1< Vector > Vectors;


public:

	LK_BallEnergy( EnergyMethodOptions const & options );


	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	LK_BallEnergy( LK_BallEnergy const & src );

	///
	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & rotamer_set
	) const;

	virtual
	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid ) const;

	virtual
 	void
	setup_for_derivatives(
		pose::Pose & pose,
		ScoreFunction const & scfxn
	) const;

	/// helper function for outside use
	Real
	calculate_lk_desolvation_of_single_atom_by_residue(
																										 Size const atom1,
																										 conformation::Residue const & rsd1,
																										 conformation::Residue const & rsd2
																										 );
	Real
	calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair(
																																	 Size const atom1,
																																	 conformation::Residue const & rsd1,
																																	 conformation::Residue const & rsd2
																																	 );
	void
	calculate_lk_ball_atom_energies(
																	Size const atom1,
																	conformation::Residue const & rsd1,
																	Vectors const & atom1_waters,
																	Size const atom2,
																	conformation::Residue const & rsd2,
																	Real & lk_desolvation_of_atom1_by_atom2,
																	Real & lk_ball_desolvation_of_atom1_by_atom2 // includes lk-fraction
																	) const;

	void
	calculate_lk_ball_atom_energies_cp(
																		 Size const atom1,
																		 conformation::Residue const & rsd1,
																		 Vectors const & atom1_waters,
																		 Size const atom2,
																		 conformation::Residue const & rsd2,
																		 etable::count_pair::CPCrossoverBehavior const & cp_crossover,
																		 Real & lk_desolvation_of_atom1_by_atom2,
																		 Real & lk_ball_desolvation_of_atom1_by_atom2 // includes lk-fraction
																		 ) const;

	// helper
	Real
	get_lk_fractional_contribution_for_single_water(
																									Vector const & atom2_xyz,
																									Size const atom2_type,
																									Vector const & atom1_water
																									) const;

	void
	eval_desolvation_derivs_no_count_pair(
																				Real const d2,
																				Size const atom1,
																				conformation::Residue const & rsd1,
																				Size const atom2,
																				conformation::Residue const & rsd2,
																				Real & atom1_lk_desolvation_by_atom2_deriv,
																				Real & atom2_lk_desolvation_by_atom1_deriv
																				);

// 	/// called during gradient-based minimization inside dfunc
// 	/**
// 		 F1 and F2 are not zeroed -- contributions from this atom are
// 		 just summed in
// 	**/
// 	virtual
// 	void
// 	eval_atom_derivative(
// 		id::AtomID const & id,
// 		pose::Pose const & pose,
// 		kinematics::DomainMap const & domain_map,
// 		ScoreFunction const & sfxn,
// 		EnergyMap const & weights,
// 		Vector & F1,
// 		Vector & F2
// 	) const;

	virtual
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
	) const;



	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	void
	residue_pair_energy(
											conformation::Residue const & rsd1,
											LKB_ResidueInfo const & rsd1_info,
											conformation::Residue const & rsd2,
											LKB_ResidueInfo const & rsd2_info,
											EnergyMap & emap
											) const;


	void
	accumulate_single_atom_contributions(
																			 Size const atom1,
																			 Size const atom1_type_index,
																			 Vectors const & atom1_waters,
																			 utility::vector1< Real > const & atom1_wts,
																			 conformation::Residue const & rsd1,
																			 Size const atom2_type_index,
																			 Vector const & atom2_xyz,
																			 Real const lk_desolvation_of_atom1_by_atom2,
																			 EnergyMap & emap
																			 ) const;

	/* Undefined, commenting out to fix PyRosetta build  void
	get_scorefxn_weights_for_derivatives(
																			 Size const atom1,
																			 conformation::Residue const & rsd1,
																			 bool const atom1_has_waters,
																			 Size const atom2,
																			 conformation::Residue const & rsd2,
																			 EnergyMap const & weights,
																			 Real & unoriented_weight,
																			 Real & oriented_weight
																			 ) const;
																			 */

	void
	setup_for_minimizing_for_residue(
																	 conformation::Residue const & rsd,
																	 pose::Pose const & pose,
																	 ScoreFunction const & scorefxn,
																	 kinematics::MinimizerMapBase const & min_map,
																	 ResSingleMinimizationData & resdata
																	 ) const;

	void
	setup_for_minimizing_for_residue_pair(
																				conformation::Residue const & rsd1,
																				conformation::Residue const & rsd2,
																				pose::Pose const &,
																				ScoreFunction const & scorefxn,
																				kinematics::MinimizerMapBase const & min_map,
																				ResSingleMinimizationData const & res1data,
																				ResSingleMinimizationData const & res2data,
																				ResPairMinimizationData & pairdata
																				) const;

	bool
	use_extended_residue_pair_energy_interface() const;

	void
	residue_pair_energy_ext(
													conformation::Residue const & rsd1,
													conformation::Residue const & rsd2,
													ResPairMinimizationData const & pairdata,
													pose::Pose const &,// pose,
													ScoreFunction const &,
													EnergyMap & emap
													) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	bool
	requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & ) const;

	void
	setup_for_scoring_for_residue(
																conformation::Residue const & rsd,
																pose::Pose const &,// pose,
																ScoreFunction const & sfxn,
																ResSingleMinimizationData & resdata
																) const;

	bool
	requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const &  ) const;

	void
	setup_for_derivatives_for_residue(
																		conformation::Residue const & rsd,
																		pose::Pose const & pose,
																		ScoreFunction const & sfxn,
																		ResSingleMinimizationData & min_data
																		) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	virtual
	Distance
	atomic_interaction_cutoff() const;


	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;



/////////////////////////////////////////////////////////////////////////////
// private methods
//private:
/////////////////////////////////////////////////////////////////////////////

	Real
	eval_lk_fraction( Real const d2_delta ) const;


	Real
	eval_d_lk_fraction_dr_over_r( Real const d2_delta ) const;


	Real
	get_lk_fractional_contribution(
																 Vector const & atom2_xyz,
																 Size const atom2_type_index,
																 Vectors const & atom1_waters,
																 Size & closest_water,
																 Real & closest_water_dis2
																 ) const;
	Real
	get_lk_fractional_contribution(
																 Vector const & atom2_xyz,
																 Size const atom2_type_index,
																 Vectors const & atom1_waters
																 ) const;
	/// for external use
	Real
	eval_lk_ball_fraction_deriv(
															Vector const & atom2_xyz,
															Size const atom2_type_index,
															Vectors const & atom1_waters,
															bool const evaluate_deriv,
															Vector & f1,
															Vector & f2
															) const;


	/* Undefined, commenting out to fix PyRosetta build  void
	residue_pair_energy(
											conformation::Residue const & rsd1,
											utility::vector1< Vectors > const & rsd1_waters,
											conformation::Residue const & rsd2,
											utility::vector1< Vectors > const & rsd2_waters,
											pose::Pose const & pose,
											ScoreFunction const &,
											EnergyMap & emap
											) const; */

	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;



	void
	sum_deriv_contributions_for_atom_pair_one_way(
																								Size const atom1,
																								conformation::Residue const & rsd1,
																								Vectors const & atom1_waters,
																								utility::vector1< Real > const & atom1_wts,
																								Size const atom2,
																								conformation::Residue const & rsd2,
																								scoring::EnergyMap const & weights,
																								Real const weight_factor,
																								Real const d2,
																								Vector & F1,
																								Vector & F2
																								) const;

	void
	sum_deriv_contributions_for_atom_pair(
																				Size const atom1,
																				conformation::Residue const & rsd1,
																				LKB_ResidueInfo const & rsd1_info,
																				Size const atom2,
																				conformation::Residue const & rsd2,
																				LKB_ResidueInfo const & rsd2_info,
																				pose::Pose const & pose,
																				scoring::EnergyMap const & weights,
																				Real const cp_weight,
																				Vector & F1,
																				Vector & F2
																				) const;


	void
	setup_d2_bounds();

// 	void
// 	add_my_score_types();

// 	/// HACK
// 	void
// 	setup_hack();

// 	bool
// 	include_residue( conformation::Residue const & rsd ) const;




/////////////////////////////////////////////////////////////////////////////
// data
private:
/////////////////////////////////////////////////////////////////////////////

	etable::EtableCAP etable_;


	/// these guys are taken from the etable
	ObjexxFCL::FArray3D< Real > const & solv1_;
	ObjexxFCL::FArray3D< Real > const & solv2_;

	ObjexxFCL::FArray3D< Real > const & dsolv1_;

	Real const safe_max_dis2_;
	Real const etable_bins_per_A2_;

	bool const use_intra_dna_cp_crossover_4_;

	static Real const ramp_width_A2_;
	utility::vector1< Real > d2_low_;
	utility::vector1< bool > atom_type_is_charged_;

	utility::vector1< Real > lk_ball_prefactor_;
	/// HACK
	//utility::vector1< Size > positions_;
	//bool include_all_dna_;

	virtual
	core::Size version() const;

};

/// this is a  helper function for hbonds
void
apply_lk_ball_fraction_weight_for_hbonds(
																				 Size const hatm,
																				 conformation::Residue const & don_rsd,
																				 Size const aatm,
																				 conformation::Residue const & acc_rsd,
																				 Vector const & hatm_xyz,
																				 Vector const & datm_xyz,
																				 Real & unweighted_energy,
																				 bool const evaluate_derivative,
																				 hbonds::HBondDerivs & hbderivs,
																				 Real & don_lk_fraction,
																				 Real & acc_lk_fraction
																				 );

}
}
}

#endif // INCLUDED_core_scoring_methods_LK_BallEnergy_HH
