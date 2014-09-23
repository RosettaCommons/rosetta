// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EnvPairPotential.fwd.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_EnvPairPotential_fwd_hh
#define INCLUDED_core_scoring_EnvPairPotential_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class CenListInfo;
typedef utility::pointer::shared_ptr< CenListInfo > CenListInfoOP;

class EnvPairPotential;

typedef utility::pointer::shared_ptr< EnvPairPotential > EnvPairPotentialOP;
typedef utility::pointer::shared_ptr< EnvPairPotential const > EnvPairPotentialCOP;

} // ns scoring
} // ns core

#endif
