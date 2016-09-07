// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/util.cc
/// @brief  Nonmember functions for evaluating some or all energy methods on residues or residue pairs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/util.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

void
eval_scsc_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1sc_centroid,
	Vector const & r2sc_centroid,
	Real const & r1sc_radius,
	Real const & r2sc_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const scsc_d2 = r1sc_centroid.distance_squared( r2sc_centroid );
	Real const scsc_radsum = r1sc_radius + r2sc_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		Real cutoff = scsc_radsum + (*iter)->atomic_interaction_cutoff();
		if ( ! (*iter)->divides_backbone_and_sidechain_energetics() || scsc_d2 < cutoff * cutoff ) {
			(*iter)->sidechain_sidechain_energy( r1, r2, pose, sfxn, emap );
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		Real cutoff = scsc_radsum + (*iter)->atomic_interaction_cutoff();
		if ( ! (*iter)->divides_backbone_and_sidechain_energetics() || scsc_d2 < cutoff * cutoff ) {
			(*iter)->sidechain_sidechain_energy( r1, r2, pose, sfxn, emap );
		}
	}
}

void
eval_bbsc_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1bb_centroid,
	Vector const & r2sc_centroid,
	Real const & r1bb_radius,
	Real const & r2sc_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const bbsc_d2 = r1bb_centroid.distance_squared( r2sc_centroid );
	Real const bbsc_radsum = r1bb_radius + r2sc_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbsc_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbsc_d2 < cutoff * cutoff ) {
				(*iter)->backbone_sidechain_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbsc_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbsc_d2 < cutoff * cutoff ) {
				(*iter)->backbone_sidechain_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
}

void
eval_bbbb_sr2b_energies(
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	Vector const & r1bb_centroid,
	Vector const & r2bb_centroid,
	Real const & r1bb_radius,
	Real const & r2bb_radius,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	Real const bbbb_d2 = r1bb_centroid.distance_squared( r2bb_centroid );
	Real const bbbb_radsum = r1bb_radius + r2bb_radius;
	for ( auto
			iter = sfxn.ci_2b_begin(), iter_end = sfxn.ci_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbbb_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbbb_d2 < cutoff * cutoff ) {
				(*iter)->backbone_backbone_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
	for ( auto
			iter = sfxn.cd_2b_begin(), iter_end = sfxn.cd_2b_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->divides_backbone_and_sidechain_energetics() ) {
			Real cutoff = bbbb_radsum + (*iter)->atomic_interaction_cutoff();
			if ( bbbb_d2 < cutoff * cutoff ) {
				(*iter)->backbone_backbone_energy( r1, r2, pose, sfxn, emap );
			}
		}
	}
}

/// @details returns the origin if there are no backbone atoms
Vector
compute_bb_centroid(
	conformation::Residue const & res
) {
	Vector bb_centroid( 0.0 );
	Size count_n_bb( 0 );
	for ( Size ii = 1; ii <= res.type().first_sidechain_atom() - 1; ++ii ) {
		++count_n_bb;
		bb_centroid += res.xyz( ii );
	}
	if ( count_n_bb != 0 ) bb_centroid /= count_n_bb;
	return bb_centroid;
}

Real
compute_bb_radius(
	conformation::Residue const & res,
	Vector const & bb_centroid
) {
	Real bb_radius = 0;
	for ( Size ii = 1; ii <= res.type().first_sidechain_atom() - 1; ++ii ) {
		Real d2 = res.xyz( ii ).distance_squared( bb_centroid );
		if ( bb_radius < d2 ) bb_radius = d2;
	}
	return bb_radius;
}

Vector
compute_sc_centroid(
	conformation::Residue const & res
) {
	Vector centroid( 0.0 );
	Size count( 0 );
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		count += 1;
		centroid += res.xyz( ii );
	}
	if ( count == 0 ) {
		return compute_bb_centroid( res );
	} else {
		centroid /= count;
		return centroid;
	}
}

Real
compute_sc_radius(
	conformation::Residue const & res,
	Vector const & centroid
) {
	Real max_d2 = 0;
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		Real d2 = res.xyz(ii).distance_squared( centroid );
		if ( d2 > max_d2 ) max_d2 = d2;
	}
	return std::sqrt( max_d2 );
}

}
}
