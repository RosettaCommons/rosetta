// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/normalmode/NormalModeRelaxMover.fwd.hh
/// @brief   initialization for NormalMode
/// @detailed
/// @author  Hahnbeom Park

#ifndef INCLUDED_protocols_normalmode_NormalModeRelaxMover_fwd_hh
#define INCLUDED_protocols_normalmode_NormalModeRelaxMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace normalmode {

class NormalModeRelaxMover;
typedef utility::pointer::owning_ptr< NormalModeRelaxMover > NormalModeRelaxMoverOP;
typedef utility::pointer::owning_ptr< NormalModeRelaxMover const > NormalModeRelaxMoverCOP;

class CartesianNormalModeMover;
typedef utility::pointer::owning_ptr< CartesianNormalModeMover > CartesianNormalModeMoverOP;
typedef utility::pointer::owning_ptr< CartesianNormalModeMover const > CartesianNormalModeMoverCOP;

class TorsionNormalModeMover;
typedef utility::pointer::owning_ptr< TorsionNormalModeMover > TorsionNormalModeMoverOP;
typedef utility::pointer::owning_ptr< TorsionNormalModeMover const > TorsionNormalModeMoverCOP;

} // normalmode
} // protocols

#endif
