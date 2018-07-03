// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/qsar/scoring_grid/GridSet.fwd.hh
/// @brief A set of related grids
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridSet_fwd_hh
#define INCLUDED_protocols_qsar_scoring_grid_GridSet_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>


// Forward
namespace protocols {
namespace qsar {
namespace scoring_grid {

class GridSet;

typedef utility::pointer::shared_ptr< GridSet > GridSetOP;
typedef utility::pointer::shared_ptr< GridSet const > GridSetCOP;

} //protocols
} //qsar
} //scoring_grid

#endif //INCLUDED_protocols_qsar_scoring_grid_GridSet_fwd_hh
