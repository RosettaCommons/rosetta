// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Implementation of the Neighborhood class.
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include <protocols/neighbor/Neighborhood.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace neighbor {

using core::Real;

static thread_local basic::Tracer TR( "protocols.neighbor" );


/// @brief Neighborhood's constructor.
///
/// @param[in] set: pose indexes of the residues whose neighborhood is to be
/// 	computed
/// @param[in] ps: the pose
/// @param[in] ngb_fun: function determining whether two residues are neighbors
/// 	of one another
///
/// @details At the end of this function, ngb_mask[i], for i in {1,...,NPS} such
/// 	that i is not a member of set, is true iff is_ngb(R_j, N_i, ps) equals
/// 	true for at least one j in {1,...,NSET}, where	R_j indicates
/// 	ps.residue(set[j]) and N_i indicates ps.residue(i); for i in {1,...,NPS}
/// 	such that i is a member of set, ngb_mask[i] equals false.
///
Neighborhood::Neighborhood(vector1<Size> const& set, Pose const& ps,
	NGB_FUN_PTR ngb_fun) : is_ngb(ngb_fun), ngb_mask(ps.total_residue(), false) {

	/// build neighbor mask
	Size NSET = set.size();
	Size NPS = ps.total_residue();
	for(Size i=1; i<=NSET; ++i) {
		Residue const& res = ps.residue(set[i]);
		for(Size j=1; j<=NPS; ++j) {
			Residue const& ngb = ps.residue(j);
			if(is_ngb(res, ngb, ps))
				ngb_mask[j] = true;
		}
	}

	/// ignore neighbors that are part of the set
	for(Size i=1; i<=NSET; ++i) {
		ngb_mask[set[i]] = false;
	}

	// collect neighbors based on mask: the ith element of ngbs is the index of
	// the ith true element in ngb_mask, namely, the pose index of the ith neigbor
	// (i=1,...,T, where T is the number of neighbors).
	for(Size i=1; i<=NPS; ++i)
		if(ngb_mask[i])
			ngbs.push_back(i);
}


/// @brief: prints the neighbor mask as a binary string
///
/// @details: the ith output digit equals 1 if the ith element in the mask is true;
///           it equals 0 if the ith element is false (i=1,...,SIZ).
///
void Neighborhood::print_ngb_mask() const {

	Size const SIZ = ngb_mask.size();
	for(Size i=1; i<=SIZ; ++i) {
		TR << (ngb_mask[i] ? "1" : "0");
	}
}


/// @brief: returns true iff r2 is a residue whose node in the energy graph is
/// 	connected to r1's.
///
/// @param[in] r1 residue r1
/// @param[in] r2 residue r2
/// @param[in] ps pose to which both residues belong
///
bool in_nrg_graph(Residue const& r1, Residue const& r2, Pose const& ps) {

	Size const i1 = r1.seqpos();
	Size const i2 = r2.seqpos();

	core::scoring::EnergyGraph const& energy_graph( ps.energies().energy_graph() );
	for ( core::graph::Graph::EdgeListConstIter
		nit = energy_graph.get_node(i1)->const_edge_list_begin(),
		nite = energy_graph.get_node(i1)->const_edge_list_end();
		nit != nite; ++nit )
		if((*nit)->get_other_ind(i1) == i2)
			return true;

	return false;
}


/// @brief: returns true iff r2's neighbor atom falls within a cutoff distance
/// 	of r1's neighbor atom.
///
/// @param[in] r1 residue r1
/// @param[in] r2 residue r2
/// @param[in] ps pose to which both residues belong
///
/// @details The squared cutoff distance is specified within the function by
/// 	constant DCUT2.
///
bool in_ngbat_sphere(Residue const& r1, Residue const& r2, Pose const&) {

	Real const DCUT2 = 36;//64; // squared Angstrom

	numeric::xyzVector<Real> ngb1 = r1.xyz(r1.nbr_atom());
	numeric::xyzVector<Real> ngb2 = r2.xyz(r2.nbr_atom());

	return ngb1.distance_squared(ngb2) <= DCUT2;
}


} // neighbor
} // protocols
