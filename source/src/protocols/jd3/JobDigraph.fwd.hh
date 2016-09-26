// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDigraph.fwd.hh
/// @brief  directed graph base class forward declarations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobDigraph_fwd_hh
#define INCLUDED_protocols_jd3_JobDigraph_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {

class JobDirectedNode;
class JobDirectedEdge;
class JobDigraph;
class JobDigraphUpdater;

typedef utility::pointer::shared_ptr< JobDigraph > JobDigraphOP;
typedef utility::pointer::shared_ptr< JobDigraph const > JobDigraphCOP;

} // namespace jd3
} // namespace protocols

#endif
