// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/TiltMoverCreator.hh
/// @brief      Tilts a protein in the membrane (Rosetta Scripts Hook)
/// @details Tilts a span, protein or part of a pose in the membrane,
///    depending on the jump number. The tilt axis is the axis
///    perpendicular to the axis connecting the embedding centers of the
///    two partners;
///    BEWARE: CANNOT USE MEMBRANE JUMP AS JUMP NUMBER!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TiltMoverCreator_hh
#define INCLUDED_protocols_membrane_TiltMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class TiltMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TiltMoverCreator_hh
