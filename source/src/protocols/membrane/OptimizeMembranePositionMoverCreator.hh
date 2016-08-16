// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Optimizes the membrane position given the high-res score function
/// @details Optimizes the membrane position given the smooth high-res score
///    function; scans the center along the normal around the initial center
///    in 0.1A steps; scans the normal in 0.2degree steps along arches
///    over the x-axis, y-axis, xy-direction, -xy-direction; outcome is
///    deterministic
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeMembranePositionMoverCreator_hh
#define INCLUDED_protocols_membrane_OptimizeMembranePositionMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class OptimizeMembranePositionMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeMembranePositionMoverCreator_hh
