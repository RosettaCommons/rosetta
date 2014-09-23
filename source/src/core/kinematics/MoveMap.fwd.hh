// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/MoveMap.fwd.hh
/// @brief  Kinematics MoveMap forward declarations header
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_MoveMap_fwd_hh
#define INCLUDED_core_kinematics_MoveMap_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace kinematics {


// Forward
class MoveMap;
typedef utility::pointer::shared_ptr< MoveMap > MoveMapOP;
typedef utility::pointer::shared_ptr< MoveMap const > MoveMapCOP;

} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_MoveMap_FWD_HH
