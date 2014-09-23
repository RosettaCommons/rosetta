// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/copy_dofs/CopyDofs.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_copydofs_CopyDofs_FWD_HH
#define INCLUDED_core_pose_copydofs_CopyDofs_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>
#include <core/id/DOF_ID.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace pose {
namespace copydofs {

	typedef utility::vector1< std::pair< id::DOF_ID, core::Real > > CopyDofsInfo;

	class CopyDofs;
	typedef utility::pointer::shared_ptr< CopyDofs > CopyDofsOP;
	typedef utility::pointer::shared_ptr< CopyDofs const > CopyDofsCOP;

} //copydofs
} //pose
} //core

#endif
