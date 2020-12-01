// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_SugarCloseEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_HH

// Unit headers
#include <core/energy_methods/RNA_SugarCloseEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


class RNA_SugarCloseEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief ctor
	RNA_SugarCloseEnergy();

	/// @brief dtor
	~RNA_SugarCloseEnergy() override;

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		core::scoring::EnergyMap & emap
	) const;


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;


	///////////////////////////////////////////////////////////////////////////////
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	/// @brief RNA_SugarCloseEnergy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	void
	setup_sugar_ring_closure_constraints( pose::Pose & pose ) const;

	void
	add_sugar_ring_closure_constraints( conformation::Residue const & rsd, core::scoring::constraints::ConstraintSet & cst_set ) const;

	// data
private:

	// Ribose closure
	Real const scale_rna_torsion_tether_;
	Real const scale_rna_torsion_sd_;
	Distance const o4prime_c1prime_bond_length_;
	Distance const o4prime_c1prime_sd_;
	core::scoring::func::HarmonicFuncOP o4prime_c1prime_dist_harm_func_;

	Real const angle_sd_;
	Real const o4prime_c1prime_c2prime_bond_angle_;
	core::scoring::func::HarmonicFuncOP o4prime_c1prime_c2prime_angle_harm_func_;
	Real const o4prime_c1prime_first_base_bond_angle_;
	core::scoring::func::HarmonicFuncOP o4prime_c1prime_first_base_angle_harm_func_;
	Real const c4prime_o4prime_c1prime_bond_angle_;
	core::scoring::func::HarmonicFuncOP c4prime_o4prime_c1prime_angle_harm_func_;

	//phenix-based constraint
	bool const use_phenix_sugar_close_;
	Distance const o4prime_c1prime_bond_north_;
	Distance const o4prime_c1prime_bond_south_;
	Distance const bond_sd_;
	Real const o4prime_c1prime_c2prime_angle_north_;
	Real const o4prime_c1prime_c2prime_angle_south_;
	Real const o4prime_c1prime_n1_9_angle_north_;
	Real const o4prime_c1prime_n1_9_angle_south_;
	Real const c4prime_o4prime_c1prime_angle_north_;
	Real const c4prime_o4prime_c1prime_angle_south_;
	Real const angle_sd1_, angle_sd2_;
	core::scoring::func::FuncOP fade_delta_north_, fade_delta_south_;

	chemical::rna::RNA_FittedTorsionInfoOP rna_fitted_torsion_info_; // currently just used for delta cutoff
	mutable core::scoring::constraints::ConstraintSetOP rna_sugar_close_constraints_;

	core::Size version() const override;


};

} //scoring
} //core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
