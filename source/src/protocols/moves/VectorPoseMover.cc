// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/VectorPoseMover.cc
/// @brief Designates a mover that can be passed multiple poses by the MSDJobDistributor
/// @brief Any movers deriving from this subclass can then act on all of the input poses simultaneously
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/moves/VectorPoseMover.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace moves {

// Default constructor does nothing
VectorPoseMover::VectorPoseMover() :
	Mover ()
{}

VectorPoseMover::VectorPoseMover( std::string const& name ) :
	Mover( name )
{}

VectorPoseMover::VectorPoseMover( VectorPoseMover const & other ) :
	Mover( other )
{}

VectorPoseMover::~VectorPoseMover() {}


std::string VectorPoseMover::get_name() const {
	return "VectorPoseMover";
}

void VectorPoseMover::set_poses( utility::vector1< core::pose::PoseOP > const & poses ) {
	poses_ = poses;
}

} //moves
} //protocols
