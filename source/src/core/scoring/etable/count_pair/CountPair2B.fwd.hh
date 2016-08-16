// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPair2B.fwd.hh
/// @brief  Count pair for residue pairs separated by another residue.
///         If you need to weight interactions between atoms separated
///         by four bonds, you could find yourself in this situation.
///         For instance, C in residue 1 and N in residue 3 are separated
///         by four bonds, and the other CountPair options don't
///         handle this correctly.
/// @author Jim Havranek (havranek@genetics.wustl.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPair2B_fwd_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPair2B_fwd_hh

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

template < class CrossoverBehavior >
class CountPair2B;

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
