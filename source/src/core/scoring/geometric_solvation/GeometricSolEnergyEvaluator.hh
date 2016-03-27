// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh
/// @brief  Geometric solvation energy.
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergyEvaluator_HH
#define INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergyEvaluator_HH

#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

#include <core/types.hh>

// Package headers
#include <core/scoring/hbonds/HBEvalTuple.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

namespace core {
namespace scoring {
namespace geometric_solvation {

class GeometricSolEnergyEvaluator: public utility::pointer::ReferenceCount {

public:

	//constructor
	GeometricSolEnergyEvaluator( methods::EnergyMethodOptions const & opts );

	GeometricSolEnergyEvaluator( GeometricSolEnergyEvaluator const & src );

	//destructor
	~GeometricSolEnergyEvaluator();

public:

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	Real
	eval_atom_energy(
		id::AtomID const & atom_id,
		pose::Pose const & pose
	) const;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & ,
		EnergyMap & emap
	) const;

	Distance
	atomic_interaction_cutoff() const;

	Real
	geometric_sol_one_way_bb_bb(
		conformation::Residue const & polar_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	Real
	geometric_sol_one_way_sc(conformation::Residue const & polar_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	Real
	res_res_geometric_sol_one_way(conformation::Residue const & polar_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	bool
	atom_is_heavy( conformation::Residue const & rsd, Size const atm ) const;


	void
	get_atom_atom_geometric_solvation_for_donor(
		Size const & don_h_atm,
		conformation::Residue const & don_rsd,
		Size const & occ_atm,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose,
		Real & energy,
		bool const update_deriv = false,
		hbonds::HBondDerivs & deriv = core::scoring::hbonds::DUMMY_DERIVS,
		hbonds::HBEvalTuple & hbe = core::scoring::hbonds::DUMMY_HBE
	) const;

	void
	get_atom_atom_geometric_solvation_for_acceptor(
		Size const & acc_atm,
		conformation::Residue const & acc_rsd,
		Size const & occ_atm,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose,
		Real & energy,
		bool const update_deriv = false,
		hbonds::HBondDerivs & deriv = core::scoring::hbonds::DUMMY_DERIVS,
		hbonds::HBEvalTuple & hbe = core::scoring::hbonds::DUMMY_HBE
	) const;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real const & geom_sol_intra_weight,
		utility::vector1< DerivVectorPair > & atom_derivs,
		bool const just_rna = false /*legacy*/
	) const;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		Real const geom_sol_weight,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	Real
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose
	) const;

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool res_moving_wrt_eachother
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const res1,
		Size const res2,
		pose::Pose const & pose
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const & res
	) const;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData & pair_data
	) const;

	Real
	precalculate_bb_bb_energy_for_design(
		pose::Pose const & pose
	) const;

private:

	bool
	check_path_distance( conformation::Residue const & rsd1, conformation::Residue const & rsd2, Size const & atm1, Size const & atm2 ) const;

	inline
	Real
	donorRes_occludingRes_geometric_sol_one_way(
		conformation::Residue const & don_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	inline
	Real
	acceptorRes_occludingRes_geometric_sol_one_way(
		conformation::Residue const & acc_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	inline
	Real
	donorRes_occludingRes_geometric_sol_one_way_bb_bb(
		conformation::Residue const & don_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	inline
	Real
	acceptorRes_occludingRes_geometric_sol_one_way_bb_bb(
		conformation::Residue const & acc_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	inline
	Real
	donorRes_occludingRes_geometric_sol_one_way_sc(
		conformation::Residue const & don_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	inline
	Real
	acceptorRes_occludingRes_geometric_sol_one_way_sc(
		conformation::Residue const & acc_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

	///////////////////////////////////////////////////////////


	inline
	Real
	occluded_water_hbond_penalty(
		bool const & is_donor,
		hbonds::HBEvalTuple const & hbond_eval_type,
		Vector const & polar_atm_xyz,
		Vector const & base_atm_xyz,
		Vector const & base2_atm_xyz,
		Vector const & occluding_atm_xyz,
		Size const & polar_nb,
		Size const & occ_nb,
		bool const update_deriv = false,
		hbonds::HBondDerivs & deriv = core::scoring::hbonds::DUMMY_DERIVS
	) const;

	inline
	void
	set_water_base_atm(
		Vector const & base_v,
		Vector const & atom_v,
		Vector const & water_v,
		Vector & water_base_v,
		Real const & xH /*cos(theta)*/,
		Distance const & bond_length ) const;

	inline
	Real
	get_water_cos(
		Vector const & base_atm_xyz,
		Vector const & polar_atm_xyz,
		Vector const & occluding_atm_xyz ) const;

	Real
	donorRes_occludingRes_geometric_sol_intra(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		bool const just_RNA = false /*legacy*/ ) const;

	Real
	acceptorRes_occludingRes_geometric_sol_intra(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		bool const just_RNA = false /*legacy*/ ) const;

	Vector
	get_acceptor_base_atm_xyz( conformation::Residue const & acc_rsd, Size const & acc_atm,
		hbonds::HBEvalTuple const & hbt ) const;

private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	methods::EnergyMethodOptions const & options_;

	// no Hbonds longer than sqrt of this (the square)
	hbonds::HBondDatabaseCOP hb_database_;
	Real const dist_cut2_;
	Real const geometric_sol_scale_;
	int const interres_path_distance_cutoff_;
	int const intrares_path_distance_cutoff_;
	mutable int path_distance_;
	bool const verbose_;

};

} //geometric_solvation
} //scoring
} //core

#endif
