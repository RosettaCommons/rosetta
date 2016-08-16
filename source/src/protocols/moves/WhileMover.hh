// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_WhileMover_hh
#define INCLUDED_protocols_moves_WhileMover_hh

// Unit headers
#include <protocols/moves/WhileMover.fwd.hh>

// Unit Headers
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace moves {

class WhileMover : public Mover {
public:
	// default constructor (nmoves=1)
	WhileMover();

	WhileMover(
		MoverOP mover_in,
		core::Size nmoves_in,
		PoseConditionOP condition
	);

	~WhileMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	MoverOP mover_;
	core::Size nmoves_;
	PoseConditionOP p_cond_;
};

class PoseCondition : public utility::pointer::ReferenceCount {
public:
	PoseCondition() {};
	virtual ~PoseCondition();
	virtual bool operator() ( core::pose::Pose const& ) = 0;
};

} // moves
} // protocols


#endif
