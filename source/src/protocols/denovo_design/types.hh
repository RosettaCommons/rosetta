// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <set>
#include <vector>

#ifdef PYROSETTA
#include <core/kinematics/FoldTree.hh>
#endif

namespace protocols {
namespace denovo_design {

static char const PARENT_DELIMETER = '.';

typedef std::string SegmentName;
typedef std::vector< SegmentName > SegmentNames;
typedef std::list< SegmentName > SegmentNameList;
typedef std::set< SegmentName > SegmentNameSet;

typedef std::pair< core::Size, core::Size > SizePair;

} // denovo_design
} // protocols

#endif
