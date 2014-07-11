// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/ppo_torsion_bin.fwd.hh
/// @brief  An enumeration to represent a binning of the phi/psi/omega torsions
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_conformation_ppo_torsion_bin_FWD_HH
#define INCLUDED_core_conformation_ppo_torsion_bin_FWD_HH

// Utility headers
#include <utility/vector0.hh>

namespace core {
namespace conformation {

/// @brief Enumeration representing the different phi/psi/omega torsion bins
enum ppo_torsion_bin {
	ppo_torbin_A = 1, // trans omega: A-G
	ppo_torbin_B,
	ppo_torbin_E,
	ppo_torbin_G,
	ppo_torbin_a, // cis omega: a-g
	ppo_torbin_b,
	ppo_torbin_e,
	ppo_torbin_g,
	ppo_torbin_X, // wild card -- "all regions" or "none of the above"
	ppo_torbin_U, // "unassigned"
	n_ppo_torsion_bins = ppo_torbin_U // keep this element last
};

typedef utility::vector0< ppo_torsion_bin > torsion_bin_string;

} // conformation
} // core



#endif
