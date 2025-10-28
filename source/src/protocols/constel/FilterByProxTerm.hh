// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/FilterByProxTerm.hh
/// @brief a class to filter out constellations based on proximity constraints
///  on the N- and C-termini of their chains.
/// @author Andrea Bazzoli

#ifndef INCLUDED_FilterByProxTerm_hh
#define INCLUDED_FilterByProxTerm_hh

#include <protocols/constel/ChainTerm.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <map>

namespace protocols {
namespace constel {

class FilterByProxTerm {

private:

	/// @brief N- and C-terminal residues for each chain in the pose
	static utility::vector1<ChainTerm> chains_;

	/// @brief number of chains in the pose
	static core::Size nchains_;

	/// @brief squared value of the maximum distance for a constellation to be
	///  considered proximal to the termini
	static core::Real max_ct_dist2_;

	/// @brief number of residues forming either terminus in a chain
	static core::Size nres_;

	/// @brief chains with termini proximal to one another
	static std::map<std::string, bool> proxnc_;

public:

	/// @brief filter initialization
	static void init(core::pose::Pose const &ps, core::Real max_ct_dist, core::Real max_tt_dist, core::Size nres);

	/// @brief tells whether a constellation satisfies the filter
	static bool is_satisfied(core::pose::Pose const &ps, utility::vector1<core::Size> const &cnl);
};

/// @brief tells whether a chain has proximal termini
bool has_prox_termini(core::pose::Pose const &ps, ChainTerm const &chain, core::Size NRES,
	core::Real DMAX2);

} // constel
} // protocols

#endif
