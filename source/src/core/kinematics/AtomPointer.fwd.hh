// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/AtomPointer.fwd.hh
/// @brief  kinematics::AtomPointer2D forward declarations header
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_kinematics_AtomPointer_fwd_hh
#define INCLUDED_core_kinematics_AtomPointer_fwd_hh


// Package headers
#include <core/id/AtomID_Map.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace kinematics {


// Forward
typedef    id::AtomID_Map< tree::AtomOP > AtomPointer2D;
typedef  utility::vector1< tree::AtomOP > AtomPointer1D; // should be the same as AtomID_Map::AtomMap


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_AtomPointer_FWD_HH
