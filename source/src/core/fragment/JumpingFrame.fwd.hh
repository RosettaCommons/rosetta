// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/fragment.fwd.hh
/// @brief  forward declarations for fragment stuff
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


#ifndef INCLUDED_core_fragment_JumpingFrame_fwd_hh
#define INCLUDED_core_fragment_JumpingFrame_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>


namespace core {
namespace fragment {

// Forward
class JumpingFrame;

typedef utility::pointer::shared_ptr< JumpingFrame > JumpingFrameOP;
typedef utility::pointer::shared_ptr< JumpingFrame const > JumpingFrameCOP;

// Forward
class NonContinuousFrame;

typedef utility::pointer::shared_ptr< NonContinuousFrame > NonContinuousFrameOP;
typedef utility::pointer::shared_ptr< NonContinuousFrame const > NonContinuousFrameCOP;


} // namespace fragment
} // namespace core

#endif
