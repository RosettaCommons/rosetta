// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Lin Jiang


#ifndef INCLUDED_protocols_moves_symmetry_SetupForFibrilMover_hh

#define INCLUDED_protocols_moves_symmetry_SetupForFibrilMover_hh

// Unit headers
#include <protocols/moves/symmetry/SetupForFibrilMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/symmetry/fibril_util.hh>

#include <core/kinematics/Jump.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace moves {
namespace symmetry {
///////////////////////////////////////////////////////////////////////////////
class SetupForFibrilMover : public Mover
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

} // symmetry
} // moves
} // rosetta
#endif
