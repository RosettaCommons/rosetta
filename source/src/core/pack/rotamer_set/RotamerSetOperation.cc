// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSetOperation.cc
/// @brief  rotamer set operation implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


// Unit Headers
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <utility/vector1.hh>


#ifdef WIN32
#include <core/pack/task/PackerTask.hh>
#endif

namespace core {
namespace pack {
namespace rotamer_set {

RotamerOperation::RotamerOperation() :
	utility::pointer::ReferenceCount()
{}

RotamerOperation::~RotamerOperation() {}


RotamerSetOperation::RotamerSetOperation() :
	utility::pointer::ReferenceCount()
{}

RotamerSetOperation::~RotamerSetOperation() {}

core::Real
RotamerSetOperation::increase_packer_residue_radius(
	pose::Pose const &, // pose,
	task::PackerTaskCOP, //the task
	core::Size //the residue index
) const
{
	return 0.0;
}

}
}
}

