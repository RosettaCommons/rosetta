// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_devel_dna_protocols_hh
#define INCLUDED_devel_dna_protocols_hh


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.fwd.hh>

//#include <protocols/moves/TrialMover.fwd.hh>


#include <utility/vector1.hh>


namespace devel {
namespace dna {


void
repack_base_pair_neighbors(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::Size const seqpos,
	bool const include_current,
	bool const repack_dna
);


void
packing_specificity_test_fast(
	core::pose::Pose const & start_pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::Size const motif_begin,
	core::Size const motif_size,
	core::Size const nloop,
	std::string const & min_type,
	core::Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const add_extra_rotamers = false,
	bool const dump_pdbs = true
);

void
packing_specificity_test_fast(
	core::pose::Pose const & start_pose,
	core::scoring::ScoreFunction const & scorefxn,
	utility::vector1< core::Size > const & motif_positions,
	core::Size const nloop,
	std::string const & min_type,
	core::Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const add_extra_rotamers = false,
	bool const dump_pdbs = true
);

/// @brief Try all possible dna basepairs at the motif positions
///evaluate their energies with a repack
///requires that the pose already have base partner info set
///
void
packing_specificity_test(
	core::pose::Pose const & start_pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::Size const motif_begin,
	core::Size const motif_size,
	std::string const & min_type,
	core::Real const min_tol,
	bool const postmin,
	std::string const & output_tag,
	bool const repack_DNA = false
);

} // namespace dna
} // namespace devel

#endif
