// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/pack_rotamers.hh
/// @brief  pack rotamers module header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_pack_rotamers_hh
#define INCLUDED_core_pack_pack_rotamers_hh

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>

// Project Headers
#include <core/types.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <utility/vector1.hh>


// #include <vector>


namespace core {
namespace pack {


void
pack_rotamers(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task
);

void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results
);

void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results,
	utility::vector1< pose::PoseOP > & pose_list
);

void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task,
	Size const nloop
);

void
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	interaction_graph::InteractionGraphBaseOP & ig
);

// PyRosetta compatible version
interaction_graph::InteractionGraphBaseOP
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets
);

void
setup_IG_res_res_weights(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig
);

/// @brief Run simulated annealing, return the energy of the best rotamer assignment
/// found, and place the best rotamers onto the input pose.
Real
pack_rotamers_run(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig,
	utility::vector0<int> rot_to_pack = utility::vector0<int>()
);

/// @brief Run simulated annealing and return the best rotamer assignment
/// found.  This function does not modify the input pose.
void
pack_rotamers_run(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy
);


void
symmetric_pack_rotamers(
  pose::Pose & pose,
  scoring::ScoreFunction const & sfxn,
  task::PackerTaskCOP task
);

void
symmetric_pack_rotamers_setup(
  pose::Pose & pose,
  scoring::ScoreFunction const & scfxn,
  task::PackerTaskCOP task,
  rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets,
  interaction_graph::InteractionGraphBaseOP & ig
);

// PyRosetta compatible version
interaction_graph::InteractionGraphBaseOP
symmetric_pack_rotamers_setup(
  pose::Pose & pose,
  scoring::ScoreFunction const & scfxn,
  task::PackerTaskCOP task,
  rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets
);

Real
symmetric_pack_rotamers_run(
  pose::Pose & pose,
  task::PackerTaskCOP task,
  rotamer_set::symmetry::SymmetricRotamerSetsCOP rotsets,
  interaction_graph::InteractionGraphBaseOP ig,
  utility::vector0< int > rot_to_pack = utility::vector0< int >()
);


} // namespace pack
} // namespace core

#endif
