// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondSet.fwd.hh
/// @brief  Hydrogen bond set for a pose forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_hbonds_HBondSet_fwd_hh
#define INCLUDED_core_scoring_hbonds_HBondSet_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace hbonds {

class HBond;

typedef utility::pointer::shared_ptr< HBond > HBondOP;
typedef utility::pointer::shared_ptr< HBond const > HBondCOP;


class HBondSet;

typedef utility::pointer::shared_ptr< HBondSet > HBondSetOP;
typedef utility::pointer::shared_ptr< HBondSet const > HBondSetCOP;


} // namespace hbonds
} // namespace scoring
} // namespace core


#endif
