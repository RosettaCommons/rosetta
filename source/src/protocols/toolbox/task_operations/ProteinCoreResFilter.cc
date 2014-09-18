// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceOperation.cc
/// @brief  TaskOperation class that finds an interface and leaves it mobile in the PackerTask
/// @author Steven Lewis smlewi@unc.edu

// Unit Headers
#include <protocols/toolbox/task_operations/ProteinCoreResFilter.hh>
#include <protocols/toolbox/task_operations/ProteinCoreResFilterCreator.hh>
#include <protocols/simple_filters/NonSequentialNeighborsFilter.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ Headers

static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.ProteinCore" );

namespace protocols {
namespace toolbox {
namespace task_operations {

ProteinCore::ProteinCore() :
	ResFilter(),
	distance_threshold_( 8.0 ),
	neighbor_cutoff_( 10 ),
	bound_( false ),
	jump_( 1 ),
	neighbor_count_cutoff_( 6 )
	{}

void
ProteinCore::parse_tag( TagCOP tag ){
  distance_threshold_ = tag->getOption< core::Real >( "distance_threshold", 8.0 ) ;
  neighbor_cutoff_ = tag->getOption< core::Size >( "neighbor_cutoff", 10 );
  bound_ = tag->getOption< bool >( "bound", false );
  jump_ = tag->getOption< core::Size >( "jump", 1 );
	neighbor_count_cutoff_ = tag->getOption< core::Size >( "neighbor_count_cutoff", 6 );

  TR<<"jump: "<<jump_<<" distance_threshold: "<<distance_threshold_<<" neighbor_cutoff: "<<neighbor_cutoff_<<" bound: "<<bound_<<" neighbor_count_cutoff_: "<<neighbor_count_cutoff_<<std::endl;
}

bool
ProteinCore::operator() ( core::pose::Pose const & pose, core::Size index ) const{
	protocols::simple_filters::NonSequentialNeighborsFilter nsnf;

	nsnf.distance_threshold( distance_threshold_ );
	nsnf.neighbor_cutoff( neighbor_cutoff_ );
  nsnf.bound( bound_ );
  nsnf.resnum( index );
	nsnf.jump( jump_ );

  core::Size const neighbor_count( (core::Size) nsnf.compute( pose ) );
	return( neighbor_count >= neighbor_count_cutoff_ );
}

core::pack::task::operation::ResFilterOP
ProteinCoreFilterCreator::create_res_filter() const{
	return new ProteinCore;
}


} //namespace protocols
} //namespace toolbox
} //namespace task_operations

