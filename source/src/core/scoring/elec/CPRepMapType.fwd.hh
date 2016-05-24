// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CPRepMapType.fwd.hh
/// @brief  Typedefs for the cached data used by Fa_ElecEnergy, stored in the ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_scoring_elec_CPRepMapType_fwd_hh
#define INCLUDED_core_scoring_elec_CPRepMapType_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <map>
#include <string>

namespace core {
namespace scoring {
namespace elec {

typedef std::map< std::string, std::map<std::string,std::string> > CPRepMapType;

typedef utility::pointer::shared_ptr< CPRepMapType > CPRepMapTypeOP;
typedef utility::pointer::shared_ptr< CPRepMapType const > CPRepMapTypeCOP;

}
}
}

#endif
