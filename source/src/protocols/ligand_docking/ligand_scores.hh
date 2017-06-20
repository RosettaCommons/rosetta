// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief
/// @author Gordon Lemmon
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_ligand_docking_ligand_scores_hh
#define INCLUDED_protocols_ligand_docking_ligand_scores_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>
#include <utility/vector1.hh>

#include <map>

namespace protocols {
namespace ligand_docking {

std::map< std::string, core::Real >
get_interface_deltas(
	char chain,
	core::pose::Pose const & after,
	const core::scoring::ScoreFunctionOP scorefxn,
	std::string const & prefix = "",
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function = nullptr
);

/// @brief Another interesting metric -- how far does the ligand centroid move?
/// Large values indicate we're outside of the intended binding site.
std::map< std::string, core::Real >
get_ligand_travel(
	char chain,
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	std::string const & prefix = ""
);

/// @details normalizaton_function will only be used if the Grids do not have their own normalization
std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	core::Size jump_id,
	core::pose::Pose const & test_pose,
	std::string const & prefix = "",
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function = nullptr
);

/// @details normalizaton_function will only be used if the Grids do not have their own normalization
std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	char chain,
	core::pose::Pose const & test_pose,
	std::string const & prefix = "",
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function = nullptr
);

/// @details normalizaton_function will only be used if the Grids do not have their own normalization
std::map< std::string, core::Real >
get_ligand_grid_scores(
	protocols::qsar::scoring_grid::GridSet const & grid_set_prototype,
	utility::vector1< core::Size > const & residues,
	std::string const & chain_label,
	core::pose::Pose const & test_pose,
	std::string const & prefix,
	protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
);

/// @brief Calculate radius of gyration for downstream non-H atoms
/// @brief Ligands tend to bind in outstretched conformations...
std::map< std::string, core::Real >
get_radius_of_gyration(
	char chain,
	core::pose::Pose const & test_pose,
	std::string const & prefix
);

std::map< std::string, core::Real >
get_ligand_RMSDs(
	char chain,
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	std::string const & prefix = ""
);

std::map< std::string, core::Real >
get_automorphic_RMSDs(
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	core::Size test_residue_id,
	core::Size ref_residue_id,
	std::string const & prefix = ""
);

std::map< std::string, core::Real >
get_multi_residue_ligand_RMSDs(
	core::pose::Pose const & test_pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & test_residue_ids,
	utility::vector1< core::Size > const & ref_residue_ids,
	char chain = 'X',
	std::string const & prefix = ""
);


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_ligand_scores_HH
