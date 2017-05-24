// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/P_AA_ABEGO3.fwd.hh
/// @brief  Amino acid probability given an ABEGO (ramachandran bins) triplet/sequence, arrays and functions
/// @author imv@uw.edu

#ifndef INCLUDED_core_scoring_P_AA_ABEGO3_fwd_hh
#define INCLUDED_core_scoring_P_AA_ABEGO3_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class P_AA_ABEGO3;

typedef  utility::pointer::shared_ptr< P_AA_ABEGO3 > P_AA_ABEGO3_OP;
typedef  utility::pointer::shared_ptr< P_AA_ABEGO3 const > P_AA_ABEGO3_COP;

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_P_AA_FWD_HH
