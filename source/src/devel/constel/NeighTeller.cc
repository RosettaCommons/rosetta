// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief A class to determine neighboring relationships between or among residues.
/// @author jk
/// @author Andrea Bazzoli

#include <devel/constel/NeighTeller.hh>
#include <basic/Tracer.hh>
#include <core/scoring/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <iomanip>
#include <list>

namespace devel {
namespace constel {

using core::Real;
using core::Size;

static basic::Tracer TR( "devel.constel.NeighTeller" );


/// @brief NeighTeller's constructor
///
/// @param[out] ref_pose pose containing the residues whose neighbors are to be
///  identified.
///
/// @remarks For speed, a scorefxn containing only fa_atr is used.
///
NeighTeller::NeighTeller(Pose& ref_pose) :

	scorefxn(core::scoring::get_score_function()),
	fa_atr_weight(0.8), interaction_score_threshold(-0.3),
	nmap(MAXNGB, utility::vector1<bool>(MAXNGB, false)) {

	(*scorefxn)(ref_pose);
}


/// @brief Tells whether a probe residue is a neighbor of a target residue.
///
/// @param[in] tgt the target residue.
/// @param[in] prb the probe residue.
/// @ref_pose[in] pose both residues belong to.
///
/// @return true if prb is a neighbor of tgt; false otherwise.
///
bool NeighTeller::isneigh(core::conformation::Residue const & tgt,
	core::conformation::Residue const & prb,
	Pose const& ref_pose) {

	unweighted_emap.clear();

	core::Vector const tgt_centroid = core::scoring::compute_sc_centroid(tgt);
	core::Vector const prb_centroid = core::scoring::compute_sc_centroid(prb);
	Real const tgt_rad = core::scoring::compute_sc_radius( tgt, tgt_centroid );
	Real const prb_rad = core::scoring::compute_sc_radius( prb, prb_centroid );

	core::scoring::eval_scsc_sr2b_energies( tgt, prb, tgt_centroid, prb_centroid, tgt_rad, prb_rad, ref_pose, *scorefxn, unweighted_emap );

	Real const weighted_fa_atr = fa_atr_weight * unweighted_emap[ core::scoring::fa_atr ];

	return weighted_fa_atr < interaction_score_threshold;
}


/// @brief: tells whether a set of residues form a tree of neighbors.
///
/// @param[in] pose indexes of the residues in the set
/// @param[in] ps the pose
///
/// @return true if the set of residues form a tree of neighbors;
///  false otherwise
///
bool NeighTeller::is_neigh_tree(utility::vector1<core::Size> const& set, Pose const& ps) {

	using core::conformation::Residue;

	Size const NRES = set.size();

	// initialize neighboring table: nmap[i][j] = true iff the ith residue in the set and
	// the jth residue in the set are neighbors (i!=j, 1<=i,j<=N))
	for ( Size i=1; i<NRES; ++i ) {
		Residue const& a = ps.residue(set[i]);
		for ( Size j=i+1; j<=NRES; ++j ) {
			Residue const& b = ps.residue(set[j]);
			nmap[j][i] = nmap[i][j] = isneigh(a, b, ps);
		}
	}

	// initialize tree of neighbors and rest of residues
	utility::vector1<Size> tree_in;
	tree_in.push_back(1);

	std::list<Size> tree_out;
	for ( Size i=2; i<=NRES; ++i ) {
		tree_out.push_back(i);
	}

	// grow tree until no other residues remain
	for ( Size i=1; i<NRES; ++i ) {
		bool found = false;
		for ( Size j=1; j<=i; ++j ) {
			Size t = tree_in[j];
			for ( auto p = tree_out.begin(), end = tree_out.end(); p!=end; ++p ) {
				if ( nmap[t][*p] ) {
					tree_in.push_back(*p);
					tree_out.erase(p);
					found = true;
					break;
				}
			}
			if ( found ) break;
		}
		if ( !found ) {
			return false;
		}
	}

	return true;
}


/// @brief Prints the neighboring table for a set of residues
///
/// @param[in] nres number of residues in the set
/// @param[out] output stream
///
/// @details Line i, column j reports nmap[i][j] (1<=i,j<=nres)
///
void NeighTeller::print_nmap(Size const nres, std::ostream& os) const {

	for ( Size i=1; i<=nres; ++i ) {
		for ( Size j=1; j<=nres; ++j ) {
			os << std::setw(2) << nmap[i][j];
		}
		os << std::endl;
	}
}


/// @brief Creates the list of residues that are neighbors of a given target
///  residue.
///
/// @param[in] tgtnum residue number of the target residue in its pose.
/// @param[out] neighs boolean mask filled in by this function to represent the
///  list of neighbors. neighs[i] will be true if the ith residue in the pose
///  is a neighbor of the target; otherwise neighs[i] will be false.
/// @param[in] ps: pose that the target and its potential neighbors belong to.
///
void mk_neigh_list(core::Size const tgtnum, utility::vector1<bool>& neighs,
	Pose& ps) {

	NeighTeller nt(ps);

	core::conformation::Residue tgtres = ps.residue(tgtnum);

	for ( Size i = 1; i <= ps.size(); ++i ) {
		if ( i != tgtnum ) {
			if ( nt.isneigh(tgtres, ps.residue(i), ps) ) {
				neighs.at(i)=true;
			}
		}
	}
}

} // constel
} // devel
