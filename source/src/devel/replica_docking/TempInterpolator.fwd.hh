// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author zhe zhang

#ifndef INCLUDED_devel_replica_docking_TempInterpolator_fwd_hh
#define INCLUDED_devel_replica_docking_TempInterpolator_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace replica_docking {

// Forward
class TempInterpolatorBase;

typedef utility::pointer::shared_ptr< TempInterpolatorBase > TempInterpolatorBaseOP;
typedef utility::pointer::shared_ptr< TempInterpolatorBase const > TempInterpolatorBaseCOP;

class TempInterpolator;

typedef utility::pointer::shared_ptr< TempInterpolator > TempInterpolatorOP;
typedef utility::pointer::shared_ptr< TempInterpolator const > TempInterpolatorCOP;

class TempFixValue;

typedef utility::pointer::shared_ptr< TempFixValue > TempFixValueOP;
typedef utility::pointer::shared_ptr< TempFixValue const > TempFixValueCOP;

} // namespace replica_docking
} // namespace devel


#endif // INCLUDED_devel_replica_docking_TempInterpolator_FWD_HH
