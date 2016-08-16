// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ResidueCoordinateChangeList.fwd.hh
/// @brief  AtomTree/Conformation communication vector class forward declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_kinematics_ResidueCoordinateChangeList_fwd_hh
#define INCLUDED_core_kinematics_ResidueCoordinateChangeList_fwd_hh

#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace kinematics {


// Forward
class ResidueCoordinateChangeList;

typedef utility::pointer::shared_ptr< ResidueCoordinateChangeList > ResidueCoordinateChangeListOP;

typedef utility::vector1< Size >                  ResidueIndexList;
typedef utility::vector1< Size >::const_iterator  ResidueListIterator;


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_ResidueCoordinateChangeList_FWD_HH
