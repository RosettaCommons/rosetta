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

/// @brief a class to filter out constellations based on proximity constraints
/// 	on the N- and C-termini of their chains.
/// @author Andrea Bazzoli

#ifndef INCLUDED_FilterByProxTerm_hh
#define INCLUDED_FilterByProxTerm_hh

#include <devel/constel/ChainTerm.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <map>

using core::pose::Pose;
using utility::vector1;
using core::Size;
using core::Real;

namespace devel {
namespace constel {

class FilterByProxTerm {

private:

	/// @brief N- and C-terminal residues for each chain in the pose
	static utility::vector1<ChainTerm> chains_;

	/// @brief number of chains in the pose
	static Size nchains_;

	/// @brief squared value of the maximum distance for a constellation to be
	/// 	considered proximal	to the termini
	static Real dct_max_2_;

	/// @brief number of residues forming either terminus in a chain
	static Size nres_;

	/// @brief chains with termini proximal to one another
	static std::map<char, bool> proxnc_;

public:

	/// @brief filter initialization
	static void init(Pose const &ps, Real dct, Real dtt, Size nres);

	/// @brief tells whether a constellation satisfies the filter
	static bool sat(Pose const &ps, vector1<Size> const &cnl);
};

/// @brief tells whether a chain has proximal termini
bool has_prox_termini(Pose const &ps, ChainTerm const &chain, Size NRES,
	Real DMAX2);

} // constel
} // devel

#endif
