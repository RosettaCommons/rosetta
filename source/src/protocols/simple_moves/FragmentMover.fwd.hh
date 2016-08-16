// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/Frame.fwd.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


#ifndef INCLUDED_protocols_simple_moves_FragmentMover_fwd_hh
#define INCLUDED_protocols_simple_moves_FragmentMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {

// Forward
class FragmentMover;

typedef utility::pointer::shared_ptr< FragmentMover > FragmentMoverOP;
typedef utility::pointer::shared_ptr< FragmentMover const > FragmentMoverCOP;

// Forward
class ClassicFragmentMover;

typedef utility::pointer::shared_ptr< ClassicFragmentMover > ClassicFragmentMoverOP;
typedef utility::pointer::shared_ptr< ClassicFragmentMover const > ClassicFragmentMoverCOP;


} // namespace simple_moves
} // namespace protocols

#endif
