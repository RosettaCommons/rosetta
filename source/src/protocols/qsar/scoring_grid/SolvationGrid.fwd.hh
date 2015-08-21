// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationGrid.fwd.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_SolvationGrid_fwd_hh
#define INCLUDED_protocols_qsar_scoring_grid_SolvationGrid_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class SolvationGrid;


typedef utility::pointer::shared_ptr<SolvationGrid> SolvationGridOP;
typedef utility::pointer::shared_ptr<SolvationGrid const> SolvationGridCOP ;

}
}
}


#endif
