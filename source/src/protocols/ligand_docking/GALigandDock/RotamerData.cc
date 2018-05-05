// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandConformer.hh
///
/// @brief  Compactly represent a docked pose by storing jump + torsion dofs + pose sidechain dofs
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>
#include <core/pose/Pose.hh>
#include <numeric/Quaternion.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <map>


namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TS( "protocols.ligand_docking.GALigandDock.RotamerData" );

float ReweightableRepEnergy::score( core::Real farep_scale ) {
	return farep_scale*fa_rep_
		+ (fa_atr_wtd_+fa_sol_wtd_+fa_elec_wtd_+lk_ball_wtd_+hbond_wtd_+oneb_wtd_+penalty_wtd_) + bonus_wtd_;
}

void ReweightableRepEnergy::reset() {
	fa_rep_ = fa_atr_wtd_ = fa_sol_wtd_ = fa_elec_wtd_ = 0.0;
	lk_ball_wtd_ = hbond_wtd_ = oneb_wtd_ = penalty_wtd_ = bonus_wtd_ = 0.0;
}


RotamerPairEnergies::RotamerPairEnergies( ) {
	nrots_ = 1; // placeholder for ligand conformer
}

void RotamerPairEnergies::clear( ) {
	nrots_ = 1; // placeholder for ligand conformer
	energy_table_.clear();
	rot_counts_.clear();
}

void
RotamerPairEnergies::add_residue( core::Size nrots_i ) {
	rot_counts_.push_back(nrots_);
	nrots_ += nrots_i;
}

void
RotamerPairEnergies::finalize ( ) {
	core::Size arraySize = nrots_*(nrots_+1)/2;
	energy_table_.resize( arraySize );
	TS << "Initialized " << arraySize*sizeof(core::Real) << " byte rotamer interaction table." << std::endl;
}

ReweightableRepEnergy &
RotamerPairEnergies::energy1b( core::Size ires, core::Size irot ) {
	core::Size ii = rot_counts_[ires] + irot;
	core::Size jj = ii;
	return energy(ii,jj);
}

ReweightableRepEnergy &
RotamerPairEnergies::energy2b( core::Size ires, core::Size irot, core::Size jres, core::Size jrot ) {
	core::Size ii = rot_counts_[ires] + irot;
	core::Size jj = rot_counts_[jres] + jrot;
	return energy(std::min(ii,jj),std::max(ii,jj));
}

ReweightableRepEnergy &
RotamerPairEnergies::energyBG( core::Size ires, core::Size irot ) {
	core::Size ii = 1;
	core::Size jj = rot_counts_[ires] + irot;
	return energy(ii,jj);
}

ReweightableRepEnergy &
RotamerPairEnergies::energy( core::Size ii, core::Size jj ) {
	core::Size index = jj*(jj-1)/2 + ii;
	return energy_table_[index];
}


}
}
}
