// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergy.hh
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_lk_dome_LK_DomePackEnergy_hh
#define INCLUDED_core_pack_guidance_scoreterms_lk_dome_LK_DomePackEnergy_hh

// Unit headers
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergy.fwd.hh>
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.fwd.hh>
#include <core/scoring/lkball/LK_DomeEnergy.hh>

// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Project headers
#include <core/types.hh>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace lk_dome {


class LK_DomePackEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	LK_DomePackEnergy( core::scoring::methods::EnergyMethodOptions const &options );


	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief LK_DomePackEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief LK_DomePackEnergy is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	void finalize_total_energy( core::pose::Pose & pose, scoring::ScoreFunction const &, scoring::EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	core::Real calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		utility::vector1< core::Size > const & current_rotamer_ids,
		core::Size const substitution_position = 0
	) const override;

	void commit_considered_substitution() override;

	/// @brief Cache data from the pose in this EnergyMethod in anticipation of packing.
	///
	void set_up_residuearrayannealableenergy_for_packing( core::pose::Pose &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn) override;

	/// @brief Clear the cached data from the pose after packing.
	///
	void clean_up_residuearrayannealableenergy_after_packing( core::pose::Pose &pose ) override;

	/// @brief Disable this energy during minimization.
	void setup_for_minimizing( pose::Pose & pose, scoring::ScoreFunction const & sfxn, kinematics::MinimizerMapBase const & minmap ) const override;

	/// @brief Re-enable this energy after minimization.
	void finalize_after_minimizing( pose::Pose & pose ) const override;

private:

	core::scoring::lkball::LK_DomeEnergyCOP lk_dome_;

	mutable LK_DomeHelperOP helper_;

};

} //lk_dome
} //guidance_scoreterms
} //pack
} //core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
