// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MissingDensityToJumpMover.cc
/// @brief  Implementation of mover that inserts a jump where there is gap in the pdb. This gap corresponds to missing density.
/// @author TJ Brunette (tjbrunette@gmail.com), May 2011

#ifndef INCLUDED_protocols_moves_MissingDensityToJumpMover_hh
#define INCLUDED_protocols_moves_MissingDensityToJumpMover_hh

// Unit headers
#include <protocols/moves/MissingDensityToJumpMover.fwd.hh>

// Unit Headers
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace moves {

class MissingDensityToJumpMover : public Mover {
public:
	// default constructor (nmoves=1)
	MissingDensityToJumpMover();

	~MissingDensityToJumpMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

};

} // moves
} // protocols


#endif
