// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/copydofs/CopyDofsInfo.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pose.copydofs.CopyDofsInfo" );

namespace core {
namespace pose {
namespace copydofs {

//Constructor
CopyDofsInfo::CopyDofsInfo()
{}

//Destructor
CopyDofsInfo::~CopyDofsInfo()
{}

void
CopyDofsInfo::clear()
{
	dofs_info_.clear();
	jumps_info_.clear();
}

void
CopyDofsInfo::push_back( std::pair< id::DOF_ID, core::Real > const & dofs_info_pair)
{
	dofs_info_.push_back( dofs_info_pair );
}

void
CopyDofsInfo::push_back( std::pair< id::AtomID, core::kinematics::Jump > const & jumps_info_pair )
{
	jumps_info_.push_back( jumps_info_pair );
}

void
CopyDofsInfo::emplace_back( std::pair< id::DOF_ID, core::Real > && dofs_info_pair)
{
	dofs_info_.emplace_back( dofs_info_pair );
}

void
CopyDofsInfo::emplace_back( std::pair< id::AtomID, core::kinematics::Jump > && jumps_info_pair )
{
	jumps_info_.emplace_back( jumps_info_pair );
}

/////////////////////////////////////////////////////////////////////
// specify dof_tolerance for speed -- changing dofs (even to the same
// value) triggers pose refold which can take some time.
//
void
CopyDofsInfo::apply_dofs( pose::Pose & pose,
	core::Real const dof_tolerance /* = 1.0e-5*/ ) const
{
	for ( auto const & dof_info : dofs_info_ ) {
		if ( dof_tolerance > 0.0 ) {
			Real const dof_value_original = pose.dof( dof_info.first );
			Real const dof_value_new      = dof_info.second;
			if ( std::abs( dof_value_original - dof_value_new ) < dof_tolerance ) continue;
		}

		pose.set_dof( dof_info.first, dof_info.second );
	}

	// could do tolerance here too
	for ( auto const & jump_info : jumps_info_ ) {
		pose.set_jump( jump_info.first, jump_info.second );
	}

}

} //copydofs
} //pose
} //core
