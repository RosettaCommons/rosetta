// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/util.hh
/// @brief  Simple utilities for computing rotamer recovery between 2 poses
/// @author JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_core_pack_util_hh
#define INCLUDED_core_pack_util_hh

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>



namespace core {
namespace pack {


// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace core::pose;

/// @brief Percentage of residues which have same rotamers
core::Real residue_rotamer_recovery( Pose & pose, Pose & ref_pose, core::Real difference=10.0 );

/// @brief Percentage of rotamers recovered
core::Real rotamer_recovery( Pose & pose, Pose & ref_pose, core::Real difference=10.0 );

/// @brief Get rotamer angle differences
/// @details Outer vector is pose length, inner vector is different chi's
utility::vector1< utility::vector1< Real > > get_rotamer_angle_diffs( Pose & pose, Pose & ref_pose );




} // namespace pack
} // namespace core

#endif
