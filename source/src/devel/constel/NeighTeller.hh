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

/// @brief A class to determine neighboring relationships between or among residues.
/// @author Andrea Bazzoli

#ifndef Included_NeighTeller_HH
#define Included_NeighTeller_HH

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.fwd.hh>
#include <iostream>

namespace devel {
namespace constel {

using core::pose::Pose;
using core::conformation::Residue;

class NeighTeller {

	core::scoring::ScoreFunctionOP scorefxn;
	core::scoring::methods::EnergyMethodOptions eopts;
	core::scoring::EnergyMap unweighted_emap;
	core::Real fa_atr_weight;
	core::Real const interaction_score_threshold;

	/// @brief: boolean table telling whether two residues in a set are neighbors,
	/// 	for sets of at most	MAXNGB residues.
	static const int MAXNGB = 10;
	utility::vector1< utility::vector1<bool> > nmap;
	void print_nmap(core::Size const nres, std::ostream& os) const;

	public:
	NeighTeller(Pose& ref_pose);

	/// @brief: tells whether a probe residue is a neighbor of a target residue.
	bool isneigh(Residue const & tgt, Residue const & prb, Pose const& ref_pose);

	/// @brief: tells whether a set of residues form a tree of neighbors.
	bool is_neigh_tree(utility::vector1<core::Size> const& set, Pose const& ps);
};

/// @brief Creates the list of residues that are neighbors of a given residue.
void mk_neigh_list(core::Size const tgtnum, utility::vector1<bool>& neighs, Pose& ps);

} // constel
} // devel

#endif
