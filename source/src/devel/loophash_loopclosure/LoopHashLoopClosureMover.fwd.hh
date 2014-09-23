// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/loophash_loopclosure/LoopHashLoopClosureMover.fwd.hh
/// @brief  LoopHashLoopClosureMover class forward delcaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_devel_loophash_loopclosure_LoopHashLoopClosureMover_fwd_HH
#define INCLUDED_devel_loophash_loopclosure_LoopHashLoopClosureMover_fwd_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace devel {
namespace loophash_loopclosure {
class  LoopHashLoopClosureMover;
typedef utility::pointer::shared_ptr< LoopHashLoopClosureMover > LoopHashLoopClosureMoverOP;
typedef utility::pointer::shared_ptr< LoopHashLoopClosureMover const > LoopHashLoopClosureMoverCOP;
typedef utility::pointer::weak_ptr< LoopHashLoopClosureMover > LoopHashLoopClosureMoverAP;
typedef utility::pointer::weak_ptr< LoopHashLoopClosureMover const > LoopHashLoopClosureMoverCAP;
}
}
#endif
