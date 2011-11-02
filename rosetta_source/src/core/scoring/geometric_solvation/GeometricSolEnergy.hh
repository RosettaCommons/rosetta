// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GeometricSolEnergy.hh
/// @brief  Hydrogen bond energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergy_hh

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergy.fwd.hh>

#include <core/types.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace geometric_solvation {

///
class GeometricSolEnergy : public methods::ContextDependentTwoBodyEnergy  {
public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;
public:

	///
	GeometricSolEnergy( methods::EnergyMethodOptions const & options );

	///@brief copy c-tor
	GeometricSolEnergy( GeometricSolEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// This evaluates everything for now,
	/// but eventually may want to split this
	/// based on backbone/backbone vs. others,
	/// as is carried out in HBondEnergy.cc
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// f1 and f2 are zeroed
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	Real
	eval_atom_energy(
		id::AtomID const & atom_id,
		pose::Pose const & pose
	) const;

	///
	//	virtual
	//	void
	//	finalize_total_energy(
	//		pose::Pose & pose,
	//		ScoreFunction const &,
	//		EnergyMap & totals
	//	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	//Real
	//hydrogen_interaction_cutoff2() const;

	///@brief GeometricSolEnergy is context sensitive
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const;

private:

	inline
	Real
	res_res_geometric_sol_one_way(
		conformation::Residue const & polar_rsd,
		conformation::Residue const & occ_rsd,
		pose::Pose const & pose ) const;

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
	occluded_water_hbond_penalty(
		bool const & is_donor,
		hbonds::HBEvalType const & hbond_eval_type,
		Vector const & polar_atm_xyz,
		Vector const & base_atm_xyz,
		Vector const & occluding_atm_xyz,
		Size const & polar_nb,
		Size const & occ_nb,
		bool const update_deriv = false,
		hbonds::HBondDerivs & deriv = hbonds::DUMMY_DERIVS
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

	bool
	atom_is_donor( conformation::Residue const & rsd, Size const atm ) const;

	bool
	atom_is_donor_h( conformation::Residue const & rsd, Size const atm ) const;

	bool
	atom_is_acceptor( conformation::Residue const & rsd, Size const atm ) const;

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
		hbonds::HBondDerivs & deriv = hbonds::DUMMY_DERIVS
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
		hbonds::HBondDerivs & deriv = hbonds::DUMMY_DERIVS
	) const;

private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	methods::EnergyMethodOptionsOP options_;

	// no Hbonds longer than sqrt of this (the square)
	hbonds::HBondDatabaseCOP hb_database_;
	Real const dist_cut2_;
	Real const geometric_sol_scale_;

	bool const verbose_;
virtual
core::Size version() const;
};

} // hbonds
} // scoring
} // core

#endif

