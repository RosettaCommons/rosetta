// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/copydofs/CopyDofsInfo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_copydofs_CopyDofsInfo_HH
#define INCLUDED_core_pose_copydofs_CopyDofsInfo_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/copydofs/CopyDofsInfo.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace core {
namespace pose {
namespace copydofs {

class CopyDofsInfo: public utility::pointer::ReferenceCount {

public:

	//constructor
	CopyDofsInfo();

	//destructor
	~CopyDofsInfo();

public:

	void
	clear();

	void
	push_back( std::pair< core::id::DOF_ID, core::Real > const & dofs_info_pair );

	void
	push_back( std::pair< core::id::AtomID, core::kinematics::Jump > const & jumps_info_pair );

	void
	emplace_back( std::pair< core::id::DOF_ID, core::Real > && dofs_info_pair );

	void
	emplace_back( std::pair< core::id::AtomID, core::kinematics::Jump > && jumps_info_pair );

	void
	apply_dofs( core::pose::Pose & pose,
		core::Real const dof_tolerance = 1.0e-5 ) const;

	utility::vector1< std::pair< core::id::DOF_ID, core::Real > > const & dofs_info() const { return dofs_info_; }
	utility::vector1< std::pair< core::id::AtomID, core::kinematics::Jump > > const & jumps_info() const { return jumps_info_; }

private:

	utility::vector1< std::pair< core::id::DOF_ID, core::Real > > dofs_info_;

	// in principle, Jumps could also be stored in dofs_info via RB1, ... RB6 DOF_Types, but we don't have nice
	//  ways to go back and forth, unfortunately.
	utility::vector1< std::pair< core::id::AtomID, core::kinematics::Jump > > jumps_info_;

};

} //copydofs
} //pose
} //core

#endif
