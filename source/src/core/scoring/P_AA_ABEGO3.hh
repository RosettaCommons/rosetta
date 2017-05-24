// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/P_AA_ABEGO3.hh
/// @brief  Amino acid probability given an ABEGO (ramachandran bins) triplet/sequence, arrays and functions
/// @author imv@uw.edu

#ifndef INCLUDED_core_scoring_P_AA_ABEGO3_hh
#define INCLUDED_core_scoring_P_AA_ABEGO3_hh

// Package headers
#include <core/scoring/types.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
//#include <core/sequence/ABEGOManager.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector0.hh>

namespace core {
namespace scoring {

class P_AA_ABEGO3 : public utility::pointer::ReferenceCount
{
private:
	utility::vector0< core::Real > p_AA_ABEGO3_;

public:
	/// @brief This constructor loads the ABEGO and AA frequencies from disk and is intended
	/// to be loaded only by the ScoringManager
	P_AA_ABEGO3();

	/// @brief Read the ss-dep amino acid probability file into P_AA_ABEGO3
	void
	read_P_AA_ABEGO3();

	/// @brief Probability energies from P(aa|abego triplet)
	core::Real
	P_AA_ABEGO3_energy( char abego_previous, char abego_current, char abego_next, chemical::AA aa ) const;
};

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_P_AA_ABEGO3
