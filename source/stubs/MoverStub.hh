// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/MoverStub.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_MoverStub_HH
#define INCLUDED_protocols_moves_MoverStub_HH

// Unit Headers
#include <protocols/moves/MoverStub.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace moves {

/// @details
class MoverStub : public protocols::moves::Mover {

public:

	/// @brief
	MoverStub(
	);

	virtual ~MoverStub();

	virtual void apply( core::pose::Pose & pose );

private:

};//end MoverStub

}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_moves_MoverStub_HH
