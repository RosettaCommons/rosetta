// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddFoldUnitMover.fwd.hh
/// @author

#ifndef INCLUDED_devel_splice_AddFoldUnitMover_fwd_hh
#define INCLUDED_devel_splice_AddFoldUnitMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace splice {
class AddFoldUnitMover;
typedef utility::pointer::shared_ptr< AddFoldUnitMover > AddFoldUnitMoverOP;
typedef utility::pointer::shared_ptr< AddFoldUnitMover const > AddFoldUnitMoverCOP;

class FoldUnitUtils;
typedef utility::pointer::shared_ptr< FoldUnitUtils > FoldUnitUtilsOP;
typedef utility::pointer::shared_ptr< FoldUnitUtils const > FoldUnitUtilsCOP;

class PoseFragmentInfo;
typedef utility::pointer::shared_ptr< PoseFragmentInfo > PoseFragmentInfoOP;
typedef utility::pointer::shared_ptr< PoseFragmentInfo const > PoseFragmentInfoCOP;

class StartFreshMover;
typedef utility::pointer::shared_ptr< StartFreshMover > StartFreshMoverOP;
typedef utility::pointer::shared_ptr< StartFreshMover const > StartFreshMoverCOP;
}
}

#endif
