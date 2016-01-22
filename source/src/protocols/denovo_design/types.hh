// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/types.hh
/// @brief Useful types for Tomponent movers
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_types_hh
#define INCLUDED_protocols_denovo_design_components_types_hh

// Unit headers

// Protocol headers

// Package headers

// Core headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <list>

namespace protocols {
namespace denovo_design {

static char const PARENT_DELIMETER = '.';

typedef std::pair< core::Size, core::Size > SizePair;
typedef std::list< std::string > StringList;
typedef utility::vector1< std::string > StringVec;

struct CutAndJump {
public:
	CutAndJump() :
		cutpoint( 0 ),
		jump( 0 ) {}
	CutAndJump( core::Size const cutval, int const jumpno ) :
		cutpoint( cutval ),
		jump( jumpno ) {}
	core::Size cutpoint;
	int jump;
	core::kinematics::FoldTreeOP ft;
};
typedef utility::vector1< CutAndJump > JumpInfo;

} // denovo_design
} // protocols

#endif
