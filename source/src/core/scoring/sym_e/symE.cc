// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/symE/symE.c
/// @brief   Implementation of symmetric design bonus class
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

// Unit headers
#include <core/scoring/sym_e/symE.hh>
#include <core/scoring/sym_e/symECreator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/DenseEnergyContainer.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.symE" );

namespace core {
namespace scoring {
namespace sym_e {


/// @details This must return a fresh instance of the symE class,
/// never an instance already in use
methods::EnergyMethodOP
symECreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new symEnergy );
}

ScoreTypes
symECreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( symE_bonus );
	return sts;
}


symEnergy::symEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new symECreator ) )
{}

methods::EnergyMethodOP symEnergy::clone() const
{
	return methods::EnergyMethodOP( new symEnergy(*this) );
}

void symEnergy::setup_for_scoring(pose::Pose &pose, const ScoreFunction &) const
{
	using namespace methods;
	if ( pose.energies().long_range_container(sym_bonus_lr) == 0 ) {
		DenseEnergyContainerOP lr_container( new DenseEnergyContainer(pose.size(), symE_bonus) );
		pose.energies().set_long_range_container(sym_bonus_lr, lr_container);
	} else {
		DenseEnergyContainerOP lr_container_copied = DenseEnergyContainerOP(utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( pose.energies().nonconst_long_range_container(sym_bonus_lr) ));
		if ( lr_container_copied->size() !=pose.size() ) {
			DenseEnergyContainerOP lr_container( new DenseEnergyContainer(pose.size(), symE_bonus) );
			pose.energies().set_long_range_container(sym_bonus_lr, lr_container);
		}
	}
	pose.update_residue_neighbors();
}

void symEnergy::setup_for_packing(pose::Pose &pose, utility::vector1< bool > const &, utility::vector1< bool > const &) const
{
	pose.update_residue_neighbors();
}

void symEnergy::indicate_required_context_graphs(utility::vector1< bool > &  ) const
{}

bool symEnergy::defines_intrares_energy(const EnergyMap & ) const
{
	return false;
}

void symEnergy::eval_intrares_energy(conformation::Residue const & , pose::Pose const & , ScoreFunction const & , EnergyMap & ) const
{}

bool symEnergy::defines_residue_pair_energy(const core::pose::Pose& , platform::Size , platform::Size ) const
{
	return true;
}

methods::LongRangeEnergyType symEnergy::long_range_type() const
{
	return methods::sym_bonus_lr;
}

void symEnergy::residue_pair_energy(
	conformation::Residue const &rsd1,
	conformation::Residue const &rsd2,
	pose::Pose const &pose,
	scoring::ScoreFunction const &,
	EnergyMap & emap
) const {
	Size rsd1Pos = rsd1.seqpos();
	Size rsd2Pos = rsd2.seqpos();
	Size totalLength = pose.size();

	int symUnits = basic::options::option[ basic::options::OptionKeys::score::symE_units ]();
	core::Real symBonus = basic::options::option[ basic::options::OptionKeys::score::symE_bonus ]();

	if ( symUnits > 0 ) {
		if ( rsd1Pos % (totalLength/symUnits) ==  rsd2Pos % (totalLength/symUnits) && rsd1Pos != rsd2Pos ) {
			if ( rsd1.name() == rsd2.name() ) {
				emap[ symE_bonus ] += symBonus;
			}
		}
	}
}

core::Size
symEnergy::version() const
{
	return 1; // Initial versioning
}

} //symE
} //scoring
} //core
