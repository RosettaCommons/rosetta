// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/md/MDBase.fwd.hh
/// @brief  initialization for MD
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_md_MDBase_fwd_hh
#define INCLUDED_protocols_md_MDBase_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace devel {
namespace md {

// Forward
class MDBase;

// Types
typedef  utility::pointer::shared_ptr< MDBase >  MDBaseOP;
typedef  utility::pointer::shared_ptr< MDBase const >  MDBaseCOP;

typedef  utility::pointer::weak_ptr< MDBase >  MDBaseAP;
typedef  utility::pointer::weak_ptr< MDBase const >  MDBaseCAP;

} // namespace devel
} // namespace md

#endif
