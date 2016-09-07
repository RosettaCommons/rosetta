// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ProlineFixMover.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_moves_ProlineFixMover_hh
#define INCLUDED_protocols_simple_moves_ProlineFixMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ProlineFixMover.fwd.hh>

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace simple_moves {

class ProlineFixMover : public protocols::moves::Mover {
public:
	ProlineFixMover()  {}
	~ProlineFixMover() override = default;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
}; // ProlineFixMover

} // moves
} // protocols


#endif
