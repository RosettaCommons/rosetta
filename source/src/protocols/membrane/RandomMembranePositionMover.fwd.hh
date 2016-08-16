// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/RandomMembranePositionMover.fwd.hh
///
/// @brief      Random membrane position mover
/// @details Make random perturbations in the position of the membrane
///    as define by some center position and normal coordinate
///    (also of fixed thickness). Rotation mover rotates normal
///    from a random angle, Translation mover picks a randomly
///    translated position from the current position.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_membrane_RandomMembranePositionMover_fwd_hh
#define INCLUDED_protocols_membrane_RandomMembranePositionMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class RandomPositionRotationMover;
typedef utility::pointer::shared_ptr< RandomPositionRotationMover > RandomPositionRotationMoverOP;
typedef utility::pointer::shared_ptr< RandomPositionRotationMover const > RandomPositionRotationMoverCOP;

class RandomPositionTranslationMover;
typedef utility::pointer::shared_ptr< RandomPositionTranslationMover > RandomPositionTranslationMoverOP;
typedef utility::pointer::shared_ptr< RandomPositionTranslationMover const > RandomPositionTranslationMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_RandomMembranePositionMover_fwd_hh
