// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/sasa_scores/sasapack.hh
/// @brief  sasapack and other scores that are normalized by residue sasa
/// @author Phil Bradley

#ifndef INCLUDED_protocols_sasa_scores_sasapack_HH
#define INCLUDED_protocols_sasa_scores_sasapack_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace sasa_scores {

/// @brief Compute residue sasa values for use in deriving and assigning sasapack-like scores
void
compute_residue_sasas_for_sasa_scores(
	core::Real const probe_radius,
	core::pose::Pose const & pose,
	utility::vector1< core::Real > & rsd_sasa
);

/// @brief Compute sasapack scores for the given pose.
/// Currently only scores non-terminal, non-disulfide, protein residues.
/// The sasapack score for a residue is the difference between its SASA with a 0.5A probe
/// and the average SASA value for that residue-type in a large set of pdb structures, conditioned
/// on the SASA with a 1.4A probe.
/// The normsasa is just the difference between a residues SASA-1.4 and the average SASA-1.4 for that residue type
///
/// @note Refitting app and python code will be checked in shortly.

void
compute_sasapack_scores(
	core::pose::Pose const & pose,
	utility::vector1< core::Real > & residue_sasapack,
	utility::vector1< core::Real > & residue_normsasa,
	core::Real & average_sasapack,
	core::Real & average_normsasa
);

/// @brief Compute normalize scores for the given pose based on average energies (hence "avgE") for pdb structures.
/// Currently only scores non-terminal, non-disulfide, protein residues.
/// The "avge" score for a residue is the difference between its per-residue score and the expected per-residue
/// score for that residue type, conditioned on the residue SASA with a 1.4A probe.
/// Right now, the following scores are excluded from the avge sum since they are often very large in native structures:
///    fa_rep, fa_dun, pro_close, omega
/// as well as paa_pp for glycine, since it's just weird. Could consider refitting these
/// The normsasa is just the difference between a residues SASA-1.4 and the average SASA-1.4 for that residue type
///
/// @note Refitting app and python code will be checked in shortly.

void
compute_avge_scores(
	core::pose::Pose const & pose,
	utility::vector1< core::Real > & residue_avge,
	utility::vector1< core::Real > & residue_normsasa,
	core::Real & average_avge,
	core::Real & average_normsasa
);


} // namespace sasa_scores
} // namespace protocols


#endif //
