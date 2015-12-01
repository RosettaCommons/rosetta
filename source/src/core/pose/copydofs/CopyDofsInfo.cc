// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/copydofs/CopyDofsInfo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.copydofs.CopyDofsInfo" );

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

/////////////////////////////////////////////////////////////////////
// specify dof_tolerance for speed -- changing dofs (even to the same
// value) triggers pose refold which can take some time.
//
void
CopyDofsInfo::apply_dofs( pose::Pose & pose,
	core::Real const dof_tolerance /* = 1.0e-5*/ ) const
{

	for ( Size n = 1; n <= dofs_info_.size(); n++ ) {

		if ( dof_tolerance > 0.0 ) {
			Real const dof_value_original = pose.dof( dofs_info_[n].first );
			Real const dof_value_new      = dofs_info_[n].second;
			if ( std::abs( dof_value_original - dof_value_new ) < dof_tolerance ) continue;
		}

		pose.set_dof( dofs_info_[n].first, dofs_info_[n].second );
	}

	for ( Size n = 1; n <= jumps_info_.size(); n++ ) { // could do tolerance here too
		pose.set_jump( jumps_info_[n].first, jumps_info_[n].second );
	}

}

} //copydofs
} //pose
} //core
