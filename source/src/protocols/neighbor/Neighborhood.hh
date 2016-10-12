// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief A class to determine the neighborhood of a set of residues.
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#ifndef Incl_prot_neig_Neighborhood_hh
#define Incl_prot_neig_Neighborhood_hh

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace neighbor {

/// @brief type of functions that determine whether two residues are neighbor
///  of one another
typedef bool (*NGB_FUN_PTR)(core::conformation::Residue const& r1, core::conformation::Residue const& r2,
	core::pose::Pose const& ps);

/// @brief A class to determine the neighborhood of a set of residues within a
///  pose
class Neighborhood : public utility::pointer::ReferenceCount {

	NGB_FUN_PTR is_ngb;

	// neighbor mask: the ith element is true iff the ith residue in the pose is
	// a neighbor of the residue set (i=1,...,N, where N is the size of the
	// pose).
	utility::vector1<bool> ngb_mask;

	// set of neighbors
	utility::vector1<Size> ngbs;

	void print_ngb_mask() const;

public:
	Neighborhood(utility::vector1<Size> const& set, core::pose::Pose const& ps, NGB_FUN_PTR ngb_fun);
	utility::vector1<Size> const& get() const {return ngbs;}
};

// pairwise neighborhood functions
bool in_nrg_graph(core::conformation::Residue const& r1, core::conformation::Residue const& r2, core::pose::Pose const& ps);
bool in_ngbat_sphere(core::conformation::Residue const& r1, core::conformation::Residue const& r2, core::pose::Pose const& ps);

} // neighbor
} // protocols

#endif
