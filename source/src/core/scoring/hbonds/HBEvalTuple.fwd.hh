// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBEvalTuple.fwd.hh
/// @brief  Tuple describing data about the donor and acceptor in a single hbond
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_hbonds_HBEvalTuple_FWD_HH
#define INCLUDED_core_scoring_hbonds_HBEvalTuple_FWD_HH

#include <core/scoring/hbonds/types.hh>

namespace core {
namespace scoring {
namespace hbonds {

class HBEvalTuple;

extern HBEvalTuple DUMMY_HBE;

} // namespace hbonds
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_hbonds_HBEvalTuple_HH
