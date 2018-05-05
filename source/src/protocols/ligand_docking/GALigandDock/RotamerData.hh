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

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_RotamerData_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_RotamerData_hh

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

// two-component scores
class ReweightableRepEnergy {
public:
	ReweightableRepEnergy() {
		reset();
	}

	float score( core::Real farep_scale );

	void reset();

	float fa_rep_;
	float fa_atr_wtd_, fa_sol_wtd_, fa_elec_wtd_, lk_ball_wtd_, hbond_wtd_, oneb_wtd_, penalty_wtd_;
	float bonus_wtd_;
};

// compact representation of allowed rotamers at each position
struct PlaceableRotamer {
	core::Size resid;
	core::Size rotno;
	float prob, prob_accum;

	ReweightableRepEnergy score;
	core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo; // precomputed lkb waters
	utility::vector1< core::Real > chis;

	core::chemical::ResidueTypeCOP restype; // pointer to an instance of the residue
};

typedef utility::vector1< PlaceableRotamer > PlaceableRotamers;

// cached rotamer-pair energies
//   - dense array
class RotamerPairEnergies {
private:
	utility::vector1< ReweightableRepEnergy > energy_table_;
	utility::vector1< core::Size > rot_counts_;
	core::Size nrots_;

public:
	RotamerPairEnergies( );

	void
	clear( );

	void
	add_residue( core::Size nrots_i );

	void
	finalize ( );

	ReweightableRepEnergy &energy1b( core::Size ires, core::Size irot );

	ReweightableRepEnergy &energy2b( core::Size ires, core::Size irot, core::Size jres, core::Size jrot );

	ReweightableRepEnergy &energyBG( core::Size ires, core::Size irot );

private:
	ReweightableRepEnergy &energy( core::Size ii, core::Size jj );
};


}
}
}

#endif
