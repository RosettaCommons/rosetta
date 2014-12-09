// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SnugDock.fwd.hh
///
/// @brief forward declaration
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

#ifndef INCLUDED_protocols_antibody_SnugDock_FWD_HH
#define INCLUDED_protocols_antibody_SnugDock_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {

class SnugDock;
typedef utility::pointer::shared_ptr< SnugDock > SnugDockOP;
typedef utility::pointer::shared_ptr< SnugDock const > SnugDockCOP;

} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_SnugDock_FWD_HH
