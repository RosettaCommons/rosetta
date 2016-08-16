// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jobdist/standard_mains.hh
///
/// @brief  Main methods that implement the typical running patterns for executables.
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_jobdist_standard_mains_hh
#define INCLUDED_protocols_jobdist_standard_mains_hh

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>

#include <protocols/jobdist/Jobs.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jobdist {


/// @brief Handles loading and writing of silent files and or PDB files - superseeds the bottom two
/// functions.

void register_options_universal_main();

int universal_main(
	protocols::moves::Mover & mover,
	float thinout_factor = 0.0
);

/// @brief Reads inputs from -s/-l/-nstruct and writes (possibly gzipped) PDB files.
/// Supplied Mover is used to transform input structure into output structure.
int main_plain_mover(
	protocols::moves::Mover & mover,
	bool random_permutation = true
);

/// @brief Reads inputs from -s/-l/-nstruct and writes (possibly gzipped) PDB files.
/// Supplied Mover is used to transform input structure into output structure.
int main_plain_pdb_mover(
	protocols::moves::Mover & mover,
	core::scoring::ScoreFunctionOP scorefxn
);


/// @brief Reads inputs from -s/-l/-nstruct and writes atomtree_diff silent files.
/// Supplied Mover is used to transform input structure into output structure.
/// undefined - commenting out to make pyrosetta compile...
//int main_atomtree_diff_mover( protocols::moves::Mover & mover,  core::scoring::ScoreFunctionOP scorefxn);


/// @brief Makes BasicJob objects from command line flags -s, -l, and -nstruct.
utility::vector1< BasicJobOP > load_s_and_l();


/// @brief Helper function to safely get current output tag that's cached in Pose.
std::string get_output_tag(core::pose::Pose const & pose);

/// @brief Helper function to safely get score_map that's cached in Pose.
std::map < std::string, core::Real > get_score_map(core::pose::Pose const & pose);

} // namespace jobdist
} // namespace protocols

#endif // INCLUDED_protocols_jobdist_standard_mains_HH
