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

#ifndef INCLUDED_core_pack_prepack_pwat_rotamers_hh
#define INCLUDED_core_pack_prepack_pwat_rotamers_hh

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/annealer/SimAnnealerBase.fwd.hh>
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/vector1.fwd.hh>
#include <numeric/xyzVector.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

struct PointDwell
{
	core::Vector xyz;
	core::Real dwell;
};

typedef utility::vector1< utility::vector1< core::Vector > > SetofSets;

namespace core {
namespace pack {

// @brief get clash check atoms and neighbor atoms
void
get_clash_and_neighbors(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	utility::vector1< core::Vector > & clash_atoms,
	utility::vector1< core::Vector > & neighbor_atoms
);


// @brief return centroid from input points
Vector
centroid(
	utility::vector1< Vector > cluster
);


// @brief bottom-up clustering based on distance from cluster centroids
void
cluster_sites(
	utility::vector1< Vector > overlap_waters,
	utility::vector1< utility::vector1< Vector > > & overlap_clusters,
	Real clus_rad
);


// @brief get all backbone pwat and sidechain lkball sites
void
get_all_possible_water_sites(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP const & rotsets,
	utility::vector1< Vector > const & clash_atoms,
	utility::vector1< Vector > const & neighbor_atoms,
	utility::vector1< utility::vector1< Vector > > & lkb_wats,
	utility::vector1< utility::vector1< Vector > > & bb_wats,
	utility::vector1< Vector > & non_pack_waters,
	utility::vector1< Vector > & all_bb_pwat,
	bool exclude_exposed
);


// @brief find intersection of lkball / backbone pwat clouds
void
find_overlap(
	utility::vector1< utility::vector1< Vector > > lkb_wats,
	utility::vector1< utility::vector1< Vector > > bb_wats,
	utility::vector1< Vector > non_pack_waters,
	utility::vector1< Vector > & overlap_waters,
	bool use_average
);

/// @brief Populate rotsets for low-resoltion water packing
/// using intersection of all possible lkball sites
void
lkb_pwat_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	SetofSets & new_pwat_rotsets,
	bool exclude_exposed,
	bool use_average
);


/// @brief Run simulated annealing, return the energy of the best rotamer assignment
/// found, and place the best rotamers onto the input pose.
Real
pack_pwat_rotamers_run(
	pose::Pose & pose,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0<int> rot_to_pack,
	utility::vector1< PointDwell > & all_rot
);


/// @brief Run simulated annealing and return the best rotamer assignment
/// found.  This function does not modify the input pose.
void
pack_pwat_rotamers_run(
	pose::Pose const & pose,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	utility::vector1< PointDwell > & all_rot
);


} // namespace pack
} // namespace core

#endif
