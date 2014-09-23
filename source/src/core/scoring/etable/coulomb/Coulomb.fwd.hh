// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/coulomb/Coulomb.fwd.hh
/// @brief  forward header for Coulomb class
/// @author Matthew O'Meara


#ifndef INCLUDED_core_scoring_etable_coulomb_Coulomb_fwd_hh
#define INCLUDED_core_scoring_etable_coulomb_Coulomb_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace etable {
namespace coulomb {

class Coulomb;
typedef utility::pointer::shared_ptr< Coulomb > CoulombOP;
typedef utility::pointer::shared_ptr< Coulomb const > CoulombCOP;

} // coulomb
} // etable
} // scoring
} // core


#endif
