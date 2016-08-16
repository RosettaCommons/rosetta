// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/EnergyMap.fwd.hh
/// @brief  Vector of scores forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_EnergyMap_fwd_hh
#define INCLUDED_core_scoring_EnergyMap_fwd_hh

namespace core {
namespace scoring {

class EMapVector;
// Depricated 5/18/2010.
//class TwoBodyEMapVector;

typedef EMapVector EnergyMap;
//typedef TwoBodyEMapVector EnergyMap;

} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_EnergyMap_FWD_HH
