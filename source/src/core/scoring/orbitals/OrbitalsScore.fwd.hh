// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/orbitals/OrbitalsScore.fwd.hh
/// @brief  Created on: Jun 2, 2010
/// @author combss

#ifndef INCLUDED_core_scoring_orbitals_ORBITALSSCORE_FWD_HH
#define INCLUDED_core_scoring_orbitals_ORBITALSSCORE_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core{
namespace scoring{
namespace orbitals {
	
class OrbitalsScore;
	
typedef  utility::pointer::shared_ptr< OrbitalsScore > OrbitalsScoreOP;
	
}
}
}

#endif /* ORBITALSSCORE_FWD_HH_ */
