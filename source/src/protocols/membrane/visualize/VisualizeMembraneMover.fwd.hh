// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/membrane/VisualizeMembraneMover.hh
///
/// @brief      Visualize Membrane Planes with Virtual Residues
/// @details    Add a set of virtual residues as a third chain to the
///    membrane pose. This tool is strictly for visualization of
///    the implicit membrane and should not be present in modeling.
///    Last Modified: 6/19/14
///
/// @author  Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_fwd_hh
#define INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace visualize {

/// @brief Add membrane planes to the pose represented by
///   2 layers of MEM virtual residues
class VisualizeMembraneMover;
typedef utility::pointer::shared_ptr< VisualizeMembraneMover > VisualizeMembraneMoverOP;
typedef utility::pointer::shared_ptr< VisualizeMembraneMover const > VisualizeMembraneMoverCOP;

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_fwd_hh
