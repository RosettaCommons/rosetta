// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/znhash/ZnHash.fwd.hh
/// @brief  Forward declaration of zinc-match hash for use in optimizing zinc coordination
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

#ifndef INCLUDED_devel_znhash_ZnHash_FWD_HH
#define INCLUDED_devel_znhash_ZnHash_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace znhash {

class ZnCoord;
class ZnHash;

class ZnMatchData;
class ZnCoordinationScorer;
class ZnCoordinationConstraint;

typedef utility::pointer::shared_ptr< ZnHash > ZnHashOP;
typedef utility::pointer::shared_ptr< ZnHash const > ZnHashCOP;

typedef utility::pointer::shared_ptr< ZnCoordinationScorer > ZnCoordinationScorerOP;
typedef utility::pointer::shared_ptr< ZnCoordinationScorer const > ZnCoordinationScorerCOP;

typedef utility::pointer::shared_ptr< ZnCoordinationConstraint > ZnCoordinationConstraintOP;
typedef utility::pointer::shared_ptr< ZnCoordinationConstraint const > ZnCoordinationConstraintCOP;

}
}


#endif
