// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mover.cc
/// @brief Method code and full headers for Mover--
/// keeps heavily-included Mover.hh small and concise to maximize compiling
/// efficiency and to make the class definitions easier to read.
/// @author

// Unit Headers
#include <protocols/moves/BoolMover.hh>

// Package headers
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/MoverStatistics.hh>

// Project headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

// tracer
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

// C++ Headers

// ObjexxFCL Headers
//#include <ObjexxFCL/string.functions.hh>

namespace protocols {
namespace moves {

/// _BoolMover

_BoolMover::_BoolMover()
	: utility::pointer::ReferenceCount(),
		type_( "_BoolMoverBase" ),
		current_tag_( "NoTag" ),
		input_pose_(0),
		native_pose_(0)
{}

_BoolMover::~_BoolMover(){}

_BoolMover::_BoolMover( std::string const & type )
	: utility::pointer::ReferenceCount(),
		type_( type ),
		current_tag_( "NoTag" ),
		input_pose_(0),
		native_pose_(0)
{}

_BoolMover::_BoolMover( _BoolMover const & other )
	: utility::pointer::ReferenceCount(),
		type_( "_BoolMoverBase" ),
		current_tag_( "NoTag" ),
		input_pose_(0),
		native_pose_(0)
{
	type_ = other.type();
	current_tag_ = other.get_current_tag();
	input_pose_ = other.get_input_pose();
	native_pose_ = other.get_native_pose();
}

core::pose::PoseCOP
_BoolMover::get_input_pose() const { return input_pose_; }

void
_BoolMover::set_input_pose( core::pose::PoseCOP pose ) { input_pose_ = pose; }

core::pose::PoseCOP
_BoolMover::get_native_pose() const { return native_pose_; }

void
_BoolMover::set_native_pose( core::pose::PoseCOP pose ) { native_pose_ = pose; }

}
}
