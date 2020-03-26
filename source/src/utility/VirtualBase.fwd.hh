// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/VirtualBase.fwd.hh
/// @brief  A Base class for virtual inheritance hierarchies
/// @author Rocco Moretti (rmorettase@gmail.com)


#ifndef INCLUDED_utility_VirtualBase_fwd_hh
#define INCLUDED_utility_VirtualBase_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {

// Forward
class VirtualBase;

typedef utility::pointer::shared_ptr< VirtualBase > VirtualBaseOP;
typedef utility::pointer::shared_ptr< VirtualBase const > VirtualBaseCOP;

} // namespace utility


#endif // INCLUDED_utility_VirtualBase_fwd_hh
