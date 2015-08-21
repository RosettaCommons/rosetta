// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  stores movable_jumps info for replica docking
/// @author Zhe Zhang

#ifndef INCLUDED_protocols_docking_RigidBodyInfo_hh
#define INCLUDED_protocols_docking_RigidBodyInfo_hh

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>


namespace protocols {
namespace docking {

class RigidBodyInfo : public utility::pointer::ReferenceCount
{
public:

	RigidBodyInfo();
	void add_jump( core::Size );
	utility::vector1< core::Size > const& movable_jumps() const { return movable_jumps_; }
	virtual ~RigidBodyInfo();

private:
	utility::vector1< core::Size > movable_jumps_;
};

} // docking
} // protocols

#endif
