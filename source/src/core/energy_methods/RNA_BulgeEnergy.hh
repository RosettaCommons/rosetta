// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_BulgeEnergy.hh
/// @brief  Score function class
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_BulgeEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_BulgeEnergy_hh

// Unit headers
#include <core/energy_methods/RNA_BulgeEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


class RNA_BulgeEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief ctor
	RNA_BulgeEnergy();

	/// @brief dtor
	~RNA_BulgeEnergy() override;

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	finalize_total_energy(
		pose::Pose &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	/// @brief RNA_BulgeEnergy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	core::Size version() const override;

	// methods
private:

	bool is_RNA_bulge( conformation::Residue const & rsd ) const;

	// data
private:

	Real const bulge_bonus_;
	bool const rna_bulge_bonus_once_per_loop_;

};

} //scoring
} //core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
