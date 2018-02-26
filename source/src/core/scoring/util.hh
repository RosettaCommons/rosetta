// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/util.hh
/// @brief  Nonmember functions for evaluating some or all energy methods on residues or residue pairs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_util_HH
#define INCLUDED_core_scoring_util_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <utility/options/OptionCollection.fwd.hh>
#include <utility/vector1.hh>

#include <string>

namespace core {
namespace scoring {

/// @brief With two bounding spheres for a pair of sidechains,
/// evaluate all the sidechain/sidechain energies.  This will
/// avoid a call to EnergyMethod E's sidechain_sidechain_energiy
/// method if a) E's atomic_interaction_cutoff + r1sc_radius +
/// r2sc_radius < dist( r1sc_centroid, r2sc_centroid ) and b)
/// E returns "true" in a call to its divides_backbone_and_-
/// sidechain_energetics() method. Both context-dependent and
/// context-independent 2-body energies are evaluated in this
/// function.
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
);

/// @brief With two bounding spheres for a backbone and a sidechain,
/// evaluate all the backbone/sidechain energies.  This will
/// avoid a call to EnergyMethod E's backbone_sidechain_energiy
/// method if either a) E's atomic_interaction_cutoff + r1bb_radius +
/// r2sc_radius < dist( r1bb_centroid, r2sc_centroid ) or b)
/// E returns "false" in a call to its divides_backbone_and_-
/// sidechain_energetics() method. The reason the call is avoided if
/// "false" is returned is that, the entirety of a residue-pair-energy
/// evaluation should be returned in the sidechain_sidechain_energy
/// evaluation, if E does not implement its own versions of the bb/bb,
/// bb/sc and sc/sc energy evaluation methods. Both context-dependent and
/// context-independent 2-body energies are evaluated in this
/// function.
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
);

/// @brief With two bounding spheres for a pair of backbones,
/// evaluate all the backbone/sidechain energies.  This will
/// avoid a call to EnergyMethod E's backbone_backbone_energiy
/// method if either a) E's atomic_interaction_cutoff + r1bb_radius +
/// r2bb_radius < dist( r1bb_centroid, r2sc_centroid ) or b)
/// E returns "false" in a call to its divides_backbone_and_-
/// sidechain_energetics() method. The reason the call is avoided if
/// "false" is returned is that, the entirety of a residue-pair-energy
/// evaluation should be returned in the sidechain_sidechain_energy
/// evaluation, if E does not implement its own versions of the bb/bb,
/// bb/sc and sc/sc energy evaluation methods. Both context-dependent and
/// context-independent 2-body energies are evaluated in this
/// function.
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
);

/// @brief Compute the average coordinate of the backbone heavy atoms
/// (aka center of mass).
Vector
compute_bb_centroid(
	conformation::Residue const & r1
);

/// @brief Given a representative point for the center of the backbone,
/// compute the largest distance of all backbone heavy atoms to that point.
Real
compute_bb_radius(
	conformation::Residue const & r1,
	Vector const & r1bb_centroid
);

/// @brief Compute the average coordiante of the sidechain atoms, (aka center of mass)
/// or, if there are no side chain heavy atoms, compute the center of mass of the
/// backbone.
Vector
compute_sc_centroid(
	conformation::Residue const & r1
);

/// @brief Given a representative point for the center of the sidechain,
/// compute the largest distance of all sidechain heavy atoms to that point.
Real
compute_sc_radius(
	conformation::Residue const & r1,
	Vector const & r1sc_centroid
);

/// @brief Check if a score function is requested with incompatible option
/// flags
/// Will return true if scorefunction is "sane" and false if not.
/// If throw_exception is true, will raise an exception instead of returning false.
bool check_score_function_sanity(
	utility::options::OptionCollection const & options,
	std::string const & scorefxn_key,
	bool throw_exception = false );

}
}

#endif

