// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ThreeHelixGraph.fwd.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_ThreeHelixGraph_FWD_HH
#define INCLUDED_ThreeHelixGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace sewing {

class ThreeHelixGraph;
typedef utility::pointer::shared_ptr< ThreeHelixGraph > ThreeHelixGraphOP;
typedef utility::pointer::shared_ptr< ThreeHelixGraph const > ThreeHelixGraphCOP;

} //sewing namespace
} //devel namespace

#endif


