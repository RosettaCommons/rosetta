// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/simple_moves/UniformPositionMover.fwd.hh
///
/// @brief      Apply a uniform (deterministic move) to a given position
/// @details	Generic movers for applying a uniform rotation or translation move
///				along a jump to a given partner in the pose.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_simple_moves_UniformPositionMover_fwd_hh
#define INCLUDED_protocols_simple_moves_UniformPositionMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

/// @brief Uniform Rotation Mover
class UniformRotationMover;
typedef utility::pointer::shared_ptr< UniformRotationMover > UniformRotationMoverOP;
typedef utility::pointer::shared_ptr< UniformRotationMover const > UniformRotationMoverCOP;

/// @brief Uniform Translation Mover
class UniformTranslationMover;
typedef utility::pointer::shared_ptr< UniformTranslationMover > UniformTranslationMoverOP;
typedef utility::pointer::shared_ptr< UniformTranslationMover const > UniformTranslationMoverCOP;

} // simple_moves
} // protocols

#endif // INCLUDED_protocols_simple_moves_UniformPositionMover_fwd_hh
