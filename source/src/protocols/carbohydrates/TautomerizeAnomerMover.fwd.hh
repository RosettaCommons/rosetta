// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/carbohydrates/TautomerizeAnomerMover.fwd.hh
/// @brief   Forward declarations for TautomerizeAnomerMover.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_FWD_HH
#define INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace carbohydrates {

/// @brief  A Mover class for tautomerizing from one anomer to another at a reducing end.
class TautomerizeAnomerMover;

typedef utility::pointer::shared_ptr< TautomerizeAnomerMover > TautomerizeAnomerMoverOP;
typedef utility::pointer::shared_ptr< TautomerizeAnomerMover const > TautomerizeAnomerMoverCOP;

}  // namespace carbohydrates
}  // namespace protocols

#endif  // INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_FWD_HH
