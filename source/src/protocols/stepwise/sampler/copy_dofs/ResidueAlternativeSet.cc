// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.copy_dofs.ResidueAlternativeSet" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

//Constructor
ResidueAlternativeSet::ResidueAlternativeSet(  utility::vector1< core::pose::PoseOP > const & pose_list,
	std::map< Size, Size > const & res_map,
	Size const representative_seqpos):
	pose_list_( pose_list ),
	res_map_( res_map ),
	representative_seqpos_( representative_seqpos )
{}

//Constructor
ResidueAlternativeSet::ResidueAlternativeSet(  utility::vector1< core::pose::PoseOP > const & pose_list,
	Size const representative_seqpos):
	pose_list_( pose_list ),
	representative_seqpos_( representative_seqpos )
{
	res_map_[ representative_seqpos ] = representative_seqpos;
}

//Destructor
ResidueAlternativeSet::~ResidueAlternativeSet(){}

utility::vector1< core::pose::PoseOP >
ResidueAlternativeSet::pose_list() const
{
	return pose_list_;
}

core::pose::PoseOP
ResidueAlternativeSet::pose( Size const n ) const
{
	return pose_list_[ n ];
}

} //copy_dofs
} //sampler
} //stepwise
} //protocols
