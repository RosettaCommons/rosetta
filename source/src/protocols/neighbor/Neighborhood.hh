// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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

	using core::pose::Pose;
	using core::conformation::Residue;
	using utility::vector1;
	using core::Size;

	/// @brief type of functions that determine whether two residues are neighbor
	/// 	of one another
	typedef bool (*NGB_FUN_PTR)(Residue const& r1, Residue const& r2,
		Pose const& ps);

	/// @brief A class to determine the neighborhood of a set of residues within a
	/// 	pose
	class Neighborhood : public utility::pointer::ReferenceCount {

		NGB_FUN_PTR is_ngb;

		// neighbor mask: the ith element is true iff the ith residue in the pose is
		// a neighbor of the residue set (i=1,...,N, where N is the size of the
		// pose).
		vector1<bool> ngb_mask;

		// set of neighbors
		vector1<Size> ngbs;

		void print_ngb_mask() const;

		public:
		Neighborhood(vector1<Size> const& set, Pose const& ps, NGB_FUN_PTR ngb_fun);
		vector1<Size> const& get() const {return ngbs;}
	};

	// pairwise neighborhood functions
	bool in_nrg_graph(Residue const& r1, Residue const& r2, Pose const& ps);
	bool in_ngbat_sphere(Residue const& r1, Residue const& r2, Pose const& ps);

} // neighbor
} // protocols

#endif
