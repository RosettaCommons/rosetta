// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util
/// @brief protocols that are specific to docking low resolution
/// @detailed
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_docking_util_hh
#define INCLUDED_protocols_docking_util_hh

#include <protocols/docking/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.fwd.hh>

#include <map>

#ifdef WIN32
#include <utility/sql_database/DatabaseSessionManager.hh>
#else
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#endif




namespace protocols {
namespace docking {

/// @brief Sets up a docking fold tree based on user-specified chains for the first and second partner
void
setup_foldtree(
	core::pose::Pose & pose,
	std::string const & partner_chainID,
	DockJumps & movable_jumps );

void
setup_foldtree(
	core::pose::Pose const & pose,
	std::string const & partner_chainID,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft);


/// @brief Sets up a docking fold tree based on looking up the
/// interface in a database
void
setup_foldtree(
	core::pose::Pose & pose,
	core::Size interface_id,
	utility::sql_database::sessionOP db_session,
	DockJumps & movable_jumps );

void
setup_foldtree(
	core::pose::Pose const & pose,
	core::Size interface_id,
	utility::sql_database::sessionOP db_session,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft);


/// @brief Sets up a docking fold tree based on a map of which chains
/// are part of which partner
void
setup_foldtree(
	core::pose::Pose & pose,
	std::map< core::Size, utility::vector1< core::Size > > const & partner_to_chains,
	DockJumps & movable_jumps );

void
setup_foldtree(
	core::pose::Pose const & pose,
	std::map< core::Size, utility::vector1< core::Size > > const & partner_to_chains,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft);


/// @brief Sets up a docking fold tree based on user-specified chains for the first and second partner
void
setup_foldtree(
	core::pose::Pose & pose,
	core::Size cutpoint,
	DockJumps & movable_jumps );

void
setup_foldtree(
	core::pose::Pose const & pose,
	core::Size cutpoint,
	DockJumps & movable_jumps,
	core::kinematics::FoldTree & ft);

	
} // docking
} // protocols

#endif
