// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionData.fwd.hh
/// @brief  Defines owning pointers for BinTransitionData class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_bin_transitions_BinTransitionData_fwd_hh
#define INCLUDED_core_scoring_bin_transitions_BinTransitionData_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace bin_transitions {

class BinTransitionData; // fwd declaration
typedef utility::pointer::shared_ptr< BinTransitionData > BinTransitionDataOP;
typedef utility::pointer::shared_ptr< BinTransitionData const > BinTransitionDataCOP;
typedef utility::vector1<BinTransitionDataOP> BinTransitionDataOPs;
typedef utility::vector1<BinTransitionDataCOP> BinTransitionDataCOPs;

} //namespare bin_transitions
} //namespace scoring
} //namespace core

#endif //INCLUDED_core_scoring_bin_transitions_BinTransitionData_fwd_hh
