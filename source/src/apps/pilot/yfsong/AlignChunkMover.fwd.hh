// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rigid/RigidBodyMotionMover.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef apps_pilot_yfsong_AlignChunkMover_FWD_HH
#define apps_pilot_yfsong_AlignChunkMover_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace apps {
	namespace pilot {

class MultiTemplateAlignChunkMover;
typedef utility::pointer::owning_ptr<MultiTemplateAlignChunkMover>       MultiTemplateAlignChunkMoverOP;
typedef utility::pointer::owning_ptr<MultiTemplateAlignChunkMover const> MultiTemplateAlignChunkMoverCOP;

class CustomFragmentMover;
typedef utility::pointer::owning_ptr<CustomFragmentMover>       CustomFragmentMoverOP;
typedef utility::pointer::owning_ptr<CustomFragmentMover const> CustomFragmentMoverCOP;
		
	} // pilot
} // apps

#endif
