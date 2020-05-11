// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DEEREnergy.fwd.hh
/// @brief  Score term for data obtained with double electron-electron resonance (DEER)
/// @details This method is initiated by assigning a weight to deer_decay and providing an input
///       text file via the option -epr_deer:input_files data.txt.
///
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_methods_DEEREnergy_fwd_hh
#define INCLUDED_core_scoring_methods_DEEREnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace methods {

class DEEREnergy;

typedef utility::pointer::shared_ptr< DEEREnergy > DEEREnergyOP;
typedef utility::pointer::shared_ptr< DEEREnergy const > DEEREnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif
