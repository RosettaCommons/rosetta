// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintEnergy.hh
/// @brief Energy method that enforces the sap_constraint
/// @details sap_constraint constrains your protein to be soluble.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergy_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergy_hh

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintEnergy.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.fwd.hh>
#include <core/pack/guidance_scoreterms/sap/SapMathConstraint.fwd.hh>

// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


class SapConstraintEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	SapConstraintEnergy( core::scoring::methods::EnergyMethodOptions const &options );


	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief SapConstraintEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief SapConstraintEnergy is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	void finalize_total_energy( core::pose::Pose & pose, scoring::ScoreFunction const &, scoring::EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, utility::vector1< core::Size > const & current_rotamer_ids, core::Size const substitution_position = 0 ) const override;

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

	utility::vector1<SapConstraintHelperOP>
	get_helpers_from_pose(
		core::pose::Pose const &pose
	) const;

	utility::vector1< std::pair< utility::vector1< std::pair< Real, SapConstraintHelperCOP > >, SapMathConstraintCOP > >
	get_math_csts_from_pose(
		core::pose::Pose const & pose,
		utility::vector1<SapConstraintHelperCOP> const & helpers
	) const;

	/// @brief Is this energy disabled (e.g. for minimization)?
	mutable bool disabled_;

	mutable utility::vector1<SapConstraintHelperOP> helpers_;

	mutable utility::vector1< std::pair< utility::vector1< std::pair< Real, SapConstraintHelperCOP > >, SapMathConstraintCOP > > math_csts_;

};

} //sap
} //guidance_scoreterms
} //pack
} //core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
