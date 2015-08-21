// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Lin Jiang


#ifndef INCLUDED_protocols_fibril_SetupForFibrilMover_hh
#define INCLUDED_protocols_fibril_SetupForFibrilMover_hh

// Unit headers
#include <protocols/fibril/SetupForFibrilMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/Jump.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace fibril {
///////////////////////////////////////////////////////////////////////////////
class SetupForFibrilMover : public moves::Mover
{

public:

	// default constructor
	SetupForFibrilMover();

	//SetupForFibrilMover( std::string const & );

	~SetupForFibrilMover();

	void
	align(
		core::pose::Pose & pose,
		core::pose::Pose & monomer_pose,
		protocols::loops::Loops core,
		protocols::loops::Loops ref_core
	);

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

};

} // fibril
} // rosetta
#endif
