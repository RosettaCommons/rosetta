// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopInserter.fwd.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_loop_creation_LoopInserter_FWD_HH
#define INCLUDED_devel_loop_creation_LoopInserter_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace loop_creation {

class LoopInserter;
typedef utility::pointer::shared_ptr< LoopInserter > LoopInserterOP;
typedef utility::pointer::shared_ptr< LoopInserter const > LoopInserterCOP;

} //loop creation
} //devel

#endif


