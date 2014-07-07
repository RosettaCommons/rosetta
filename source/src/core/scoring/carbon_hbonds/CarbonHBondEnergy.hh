// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CarbonHBondEnergy.hh
/// @brief  Hydrogen bond energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_carbon_hbonds_CarbonHBondEnergy_hh
#define INCLUDED_core_scoring_carbon_hbonds_CarbonHBondEnergy_hh

#include <core/types.hh>

// Unit Headers
#include <core/scoring/carbon_hbonds/CarbonHBondEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/carbon_hbonds/CarbonHBondPotential.hh>
#include <core/scoring/MinimizationData.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>

#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace carbon_hbonds {

static core::Vector DUMMY_VECTOR;

///
class CarbonHBondEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:

	///
	CarbonHBondEnergy();

	///@brief copy c-tor
	CarbonHBondEnergy( CarbonHBondEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Evaluate the cbond energy between two residues
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	/// @brief Splitting out based on bb/bb for the OnTheFly IGs
	virtual
	void
	backbone_backbone_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	virtual
	void
	backbone_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	sidechain_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/* DEPRECATED
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
	) const;*/

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief Evaluate all atom-pair derivatives for any interactions between the two residues
	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	virtual
	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
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

	virtual
	bool
	divides_backbone_and_sidechain_energetics() const
	{
		return true;
	}

	// Real hydrogen_interaction_cutoff2() const;

	///@brief CarbonHBondEnergy is context sensitive
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const;

	bool
	use_orientation_dep_rna_ch_o_bonds(conformation::Residue const & don_rsd, conformation::Residue const & acc_rsd) const;

	/// @details Atom pair interaction energy; returns true if the interaction energy
	/// is nonzero.  The f2 vector returned when update_deriv is true is the force on
	/// the hydrogen atom; f2 for the acceptor is -1 * f2 for the hydrogen.  f1 may
	/// be computed by taking the cross product of f2 with the coordinate of the
	/// acceptor/hydrogen respectively.
	bool
	get_atom_atom_carbon_hbond_energy(
	  Size const don_h_atm,
		conformation::Residue const & don_rsd,
		Size const acc_atm,
		conformation::Residue const & acc_rsd,
		Real & energy,
		bool const update_deriv = false,
		Vector & f2 = DUMMY_VECTOR
	) const;

	Real max_dis2() const{ return max_dis2_; }

private:

	bool
	path_distance_OK(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Size const ii,
		Size const jj
	) const;

	Real
	res_res_carbon_hbond_one_way(
		conformation::Residue const & don_rsd,
		conformation::Residue const & acc_rsd,
		Real & bb_bb,
		Real & bb_sc,
		Real & sc_sc
	) const;

	Real
	bb_bb_carbon_hbond_one_way(
		conformation::Residue const & don_rsd,
		conformation::Residue const & acc_rsd
	) const;

	Real
	sc_bb_carbon_hbond_one_way(
		conformation::Residue const & don_rsd, // sidechain atoms on donor
		conformation::Residue const & acc_rsd  // backbone atoms on acceptor
	) const;

	Real
	bb_sc_carbon_hbond_one_way(
		conformation::Residue const & don_rsd, // backbone atoms on donor
		conformation::Residue const & acc_rsd  // sidechain atoms on acceptor
	) const;

	Real
	sc_sc_carbon_hbond_one_way(
		conformation::Residue const & don_rsd,
		conformation::Residue const & acc_rsd
	) const;


	void
	res_res_carbon_hbond_derivs_one_way(
		conformation::Residue const & don_rsd,
		conformation::Residue const & acc_rsd,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & don_atom_derivs,
		utility::vector1< DerivVectorPair > & acc_atom_derivs
	) const;


	bool
	atom_is_apolar_h( conformation::Residue const & rsd, Size const atm ) const;

	bool
	atom_is_acceptor( conformation::Residue const & rsd, Size const atm ) const;


private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

	// const-ref to scoring database
	carbon_hbonds::CarbonHBondPotential const & carbon_hbond_potential_;
	Real const max_dis_;
	Real const max_dis2_;

	Size const path_dist_cutoff_;

	bool const orientation_dep_rna_ch_o_bonds_;
	bool const verbose_;

	/// apl -- enmethods should not know what weights to use; they should use the weights they're given
	/// This will allow (near future) reweighting between arbitrary residue pairs
	//mutable Real wbb_bb_;
	//mutable Real wbb_sc_;
	//mutable Real wsc_sc_;

virtual
core::Size version() const;
};

} // carbon_hbonds
} // scoring
} // core

#endif

