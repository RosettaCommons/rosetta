// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSSingleGrid.fwd.hh
/// @brief   forward declaration for class PCSSingleGrid
/// @details last Modified: 05/17/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_fwd_hh
#define INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class PCSSingleGrid;

typedef utility::pointer::shared_ptr<PCSSingleGrid> PCSSingleGridOP;
typedef utility::pointer::shared_ptr<PCSSingleGrid const> PCSSingleGridCOP ;

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols

#endif //INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_fwd_hh
