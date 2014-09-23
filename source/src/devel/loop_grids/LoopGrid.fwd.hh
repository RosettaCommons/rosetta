/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/LoopGrid.fwd.hh
/// @brief  generate accessible grids for loop
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_devel_loop_grids_LoopGrid_fwd_hh
#define INCLUDED_devel_loop_grids_LoopGrid_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {

class LoopGrid;
typedef utility::pointer::shared_ptr< LoopGrid > LoopGridOP;
typedef utility::pointer::shared_ptr< LoopGrid const > LoopGridCOP;

}//namespace match
}//namespace protocols

#endif

