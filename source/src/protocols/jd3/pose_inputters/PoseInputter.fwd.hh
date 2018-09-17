// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/PoseInputter.fwd.hh
/// @brief  Forward declaration of the %PoseInputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_inputters_PoseInputter_FWD_HH
#define INCLUDED_protocols_jd3_pose_inputters_PoseInputter_FWD_HH

// Package headers
#include <protocols/jd3/pose_inputters/PoseInputSource.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

class PoseInputter;

typedef utility::pointer::shared_ptr< PoseInputter > PoseInputterOP;
typedef utility::pointer::shared_ptr< PoseInputter const > PoseInputterCOP;

typedef utility::vector1< std::pair< PoseInputSourceOP, PoseInputterOP > > PoseInputSourcesAndInputters;

} // namespace pose_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PoseInputter_FWD_HH
