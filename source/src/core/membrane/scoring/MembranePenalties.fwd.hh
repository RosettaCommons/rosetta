// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembranePenalties.fwd.hh
///
/// @brief  	Membrane Penalties
/// @detail		This class's methods evaluate scoring penalties for the membrane region
///				This will later get instantiated in MembraneSearch and MembranePtoential
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembranePenalties_fwd_hh
#define INCLUDED_core_membrane_scoring_MembranePenalties_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace scoring {

/// @brief  Class: MembranePenalties
/// @details	This class's methods evaluate scoring penalties for the membrane region
class MembranePenalties;
typedef utility::pointer::owning_ptr< MembranePenalties > MembranePenaltiesOP;
typedef utility::pointer::owning_ptr< MembranePenalties const > MembranePenaltiesCOP;

} // scoring
} // membrane
} // core

#endif // INCLUDED_core_membrane_scoring_MembranePenalties_fwd_hh