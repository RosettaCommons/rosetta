// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Zhe Zhang

#include <protocols/docking/RigidBodyInfo.hh>
#include <protocols/docking/types.hh>
//#include <iostream>


namespace protocols {
namespace docking {

RigidBodyInfo::RigidBodyInfo() {}

RigidBodyInfo::~RigidBodyInfo() {}

void RigidBodyInfo::add_jump( core::Size jump_id )
{
	movable_jumps_.push_back( jump_id );
}

} // docking
} // protocols
