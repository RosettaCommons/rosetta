// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/environment/claims/EnvLabelSelector.fwd.hh
/// @brief  The EnvLabelSelector holds a selection that another object can set inside of it.
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_claims_EnvLabelSelector_FWD_HH
#define INCLUDED_protocols_environment_claims_EnvLabelSelector_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace environment {
namespace claims {

class EnvLabelSelector;

typedef utility::pointer::shared_ptr< EnvLabelSelector > EnvLabelSelectorOP;
typedef utility::pointer::shared_ptr< EnvLabelSelector const > EnvLabelSelectorCOP;

} //namespace claims
} //namespace environment
} //namespace protocols


#endif
