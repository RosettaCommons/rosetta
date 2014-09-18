// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandDockProtocol.hh
///
/// @brief
/// @author Gordon Lemmon


#ifndef INCLUDED_protocols_ligand_docking_ligand_scores_hh
#define INCLUDED_protocols_ligand_docking_ligand_scores_hh

// AUTO-REMOVED #include <protocols/jd2/Job.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <protocols/qsar/scoring_grid/GridManager.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>

#include <core/types.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer ligand_scores_tracer( "protocols.ligand_docking.ligand_scores", basic::t_debug );

///@brief append interface_delta scores
void
append_interface_deltas(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		const core::scoring::ScoreFunctionOP scorefxn,
		std::string const & prefix
);

	
void
append_interface_deltas(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		const core::scoring::ScoreFunctionOP scorefxn,
		std::string const & prefix,
		protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
);

void
append_ligand_travel(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
);

void
append_ligand_grid_scores(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		std::string const & prefix
);

void
append_ligand_grid_scores(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & after,
		std::string const & prefix,
		protocols::qsar::scoring_grid::ScoreNormalizationOP normalization_function
);

void
append_radius_of_gyration(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		std::string const & prefix
);

void
append_ligand_RMSD(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
);

void
append_multi_residue_ligand_RMSD(
		core::Size jump_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
);

void
append_automorphic_rmsd(
		core::Size ligand_residue_id,
		protocols::jd2::JobOP job,
		core::pose::Pose const & before,
		core::pose::Pose const & after,
		std::string const & prefix
);

} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_ligand_scores_HH
