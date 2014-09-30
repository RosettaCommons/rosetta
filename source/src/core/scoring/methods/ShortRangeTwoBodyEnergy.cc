// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ShortRangeTwoBodyEnergy.cc
/// @brief  Short-Range, Two-Body Energy Method base class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {

ShortRangeTwoBodyEnergy::ShortRangeTwoBodyEnergy( EnergyMethodCreatorOP creator ) : parent( creator ) {}

ShortRangeTwoBodyEnergy::~ShortRangeTwoBodyEnergy() {}

bool
ShortRangeTwoBodyEnergy::divides_backbone_and_sidechain_energetics() const
{
	return false;
}

void
ShortRangeTwoBodyEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;

	EnergyMap emap;

	for ( Size ii = 1; ii <= set1.get_n_residue_types(); ++ii ) {
		Size const ii_offset = set1.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( set1.rotamer_ref( ii_offset ));
		Vector const & ii_coord( ii_example_rotamer.atom( ii_example_rotamer.type().nbr_atom() ).xyz());
		Real const ii_radius( ii_example_rotamer.type().nbr_radius() );

		for ( Size jj = 1; jj <= set2.get_n_residue_types(); ++jj ) {
			Size const jj_offset = set2.get_residue_type_begin( jj );
			Residue const & jj_example_rotamer( set2.rotamer_ref( jj_offset ));
			Vector const & jj_coord( jj_example_rotamer.atom( jj_example_rotamer.type().nbr_atom() ).xyz());
			Real const jj_radius( jj_example_rotamer.type().nbr_radius() );

			if ( ii_coord.distance_squared( jj_coord ) < std::pow(ii_radius+jj_radius+atomic_interaction_cutoff(), 2 )) {
				for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
					Size const kk_rot_id = ii_offset + kk - 1;
					for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
						Size const ll_rot_id = jj_offset + ll - 1;

						emap.zero( score_types() );
						residue_pair_energy( set1.rotamer_ref( kk_rot_id ), set2.rotamer_ref( ll_rot_id ), pose, sfxn, emap );
						energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy > (weights.dot( emap, score_types() ));
					}
				}
			}
		}
	}

}

void
ShortRangeTwoBodyEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero( score_types() );
		residue_pair_energy( set.rotamer_ref( ii ), residue, pose, sfxn, emap );
		energy_vector[ ii ] += static_cast< core::PackerEnergy > (weights.dot( emap, score_types() ));
	}
}

void
ShortRangeTwoBodyEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & ,	// weights
	utility::vector1< EnergyMap > & emaps
) const
{
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		emap.zero();
		residue_pair_energy( set.rotamer_ref( ii ), residue, pose, sfxn, emap );
		emaps[ ii ] += emap;
	}
}

}
}
}

