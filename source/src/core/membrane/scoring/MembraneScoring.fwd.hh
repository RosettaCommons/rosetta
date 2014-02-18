// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembraneScoring.fwd.hh
///
/// @brief  	Membrane Scoring
/// @detail		This class contains all of the independent scoring methods
///				to apply such methods to a membrane bound pose
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembraneScoring_fwd_hh
#define INCLUDED_core_membrane_scoring_MembraneScoring_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace scoring {

/// @brief  Class: Membrane Scoring
/// @detail	This class contains all of the independent scoring methods
///			to apply such methods to a membrane bound pose
class MembraneScoring;
typedef utility::pointer::owning_ptr< MembraneScoring > MembraneScoringOP;
typedef utility::pointer::owning_ptr< MembraneScoring const > MembraneScoringCOP;

} // scoring
} // membrane
} // core

#endif // INCLUDED_core_membrane_scoring_MembraneScoring_fwd_hh