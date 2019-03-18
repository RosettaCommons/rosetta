// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/VectorPoseFilter.cc
/// @brief Designates a filter that can be passed multiple poses by the VectorPoseJobDistributor
/// Any filters deriving from this subclass can then act on all of the input poses simultaneously
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/filters/VectorPoseFilter.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace filters {

// Default constructor does nothing
VectorPoseFilter::VectorPoseFilter() :
	Filter ()
{}

VectorPoseFilter::VectorPoseFilter( std::string const& name ) :
	Filter( name )
{}

VectorPoseFilter::VectorPoseFilter( VectorPoseFilter const & other ) :
	Filter( other )
{}

VectorPoseFilter::~VectorPoseFilter() {}

void VectorPoseFilter::set_poses( utility::vector1< core::pose::PoseOP > const & poses ) {
	poses_ = poses;
}

} //filters
} //protocols
