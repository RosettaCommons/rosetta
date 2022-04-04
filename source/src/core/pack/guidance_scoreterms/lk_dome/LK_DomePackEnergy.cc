// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergy.cc
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergy.hh>
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomePackEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <utility/numbers.hh>
#include <core/select/residue_selector/ResidueSelector.hh>


// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

using namespace core::scoring;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace lk_dome {

static basic::Tracer TR("core.pack.guidance_scoreterms.lk_dome.LK_DomePackEnergy");


core::scoring::methods::EnergyMethodOP
LK_DomePackEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return utility::pointer::make_shared<LK_DomePackEnergy>( options );
}

ScoreTypes
LK_DomePackEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( scoring::lk_dome_pack );
	return sts;
}

LK_DomePackEnergy::LK_DomePackEnergy ( core::scoring::methods::EnergyMethodOptions const & options ) :
	parent1( utility::pointer::make_shared<LK_DomePackEnergyCreator>() ),
	parent2( ),
	helper_()
{
	lk_dome_ = utility::pointer::make_shared<core::scoring::lkball::LK_DomeEnergy>( options );
}


core::scoring::methods::EnergyMethodOP LK_DomePackEnergy::clone() const {
	return utility::pointer::make_shared< LK_DomePackEnergy >( *this );
}

void LK_DomePackEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	return;
}

core::Size LK_DomePackEnergy::version() const
{
	return 1; // Initial versioning
}

void
LK_DomePackEnergy::finalize_total_energy(
	core::pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const {
}

core::Real
LK_DomePackEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	utility::vector1< core::Size > const & current_rotamer_ids,
	core::Size const substitution_position
) const {
	return helper_->calculate_energy( resvect, current_rotamer_ids, substitution_position );
}

void
LK_DomePackEnergy::commit_considered_substitution() {
	helper_->commit_considered_substitution();
}



/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
LK_DomePackEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose &pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets,
	core::scoring::ScoreFunction const & sfxn
) {
	TR.Debug << "Dome setup" << std::endl;

	EnergyMap const & weights = sfxn.weights();

	helper_ = utility::pointer::make_shared<LK_DomeHelper>(
		lk_dome_,
		weights[core::scoring::lk_dome],
		weights[core::scoring::lk_dome_iso],
		weights[core::scoring::lk_dome_bridge],
		weights[core::scoring::lk_dome_bridge_uncpl],
		weights[core::scoring::lk_ball_bridge2],
		weights[core::scoring::lk_ball_bridge_uncpl2]
	);

	helper_->init_with_pose( pose, rotamer_sets );

	TR.Debug << "Dome setup done" << std::endl;
}

/// @brief Clear the cached data from the pose after packing.
///
void
LK_DomePackEnergy::clean_up_residuearrayannealableenergy_after_packing(
	core::pose::Pose &
) {
	helper_ = nullptr;
}

void
LK_DomePackEnergy::setup_for_minimizing(
	pose::Pose & /*pose*/, ScoreFunction const & /*sfxn*/,
	kinematics::MinimizerMapBase const & /*minmap*/
) const {
}

void
LK_DomePackEnergy::finalize_after_minimizing( pose::Pose & /*pose*/ ) const {
}



} //lk_dome
} //guidance_scoreterms
} //pack
} //core
