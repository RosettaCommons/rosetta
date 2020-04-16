// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/pack_rotamers.hh
/// @brief  pack rotamers module header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) to allow multi-threaded interaction graph setup.

#ifndef INCLUDED_core_pack_pack_rotamers_hh
#define INCLUDED_core_pack_pack_rotamers_hh

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>

// Project Headers
#include <core/types.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {


/// @brief The entry point to calling the packer, which optimizes side-chain identities and conformations
/// (if you're designing) or conformations alone (if you're predicting structures).
/// @details Wraps the two very distinct and separate stages of rotamer packing,
/// which are factored so that they may be called asynchronously.
/// Use this wrapper as a base model for higher-level packing routines (such as pack_rotamers_loop).
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

/// @brief Run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results.
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

/// @brief Get rotamers, compute energies.
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	interaction_graph::AnnealableGraphBaseOP & ig,
	core::Size nloop = 1 //how many packing runs will use this IG? This totally optional parameter helps the factory decide which IG to use if nothing else is specified
);

/// @brief PyRosetta compatible version.
interaction_graph::AnnealableGraphBaseOP
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets
);

/// @brief Run simulated annealing, return the energy of the best rotamer assignment
/// found, and place the best rotamers onto the input pose.
Real
pack_rotamers_run(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0<int> rot_to_pack = utility::vector0<int>()
);

/// @brief Run simulated annealing and return the best rotamer assignment
/// found.  This function does not modify the input pose.
void
pack_rotamers_run(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy
);

/// @brief Provide the opportunity to clean up cached data from the pose or scorefunction after packing.
/// @details This should be called after pack_rotamers_run.  It is called from the pack_rotamers() function.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
pack_rotamers_cleanup(
	core::pose::Pose & pose,
	core::pack::interaction_graph::AnnealableGraphBaseOP annealable_graph
);

} // namespace pack
} // namespace core

#endif
