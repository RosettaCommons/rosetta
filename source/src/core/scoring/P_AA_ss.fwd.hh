// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/P_AA.hh
/// @brief  Amino acid probability arrays and functions
/// @author FD

#ifndef INCLUDED_core_scoring_P_AA_ss_fwd_hh
#define INCLUDED_core_scoring_P_AA_ss_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class P_AA_ss;

typedef  utility::pointer::shared_ptr< P_AA_ss > P_AA_ssOP;
typedef  utility::pointer::shared_ptr< P_AA_ss const > P_AA_ssCOP;

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_P_AA_FWD_HH
