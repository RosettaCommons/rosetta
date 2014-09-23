// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/SilentTrajectoryRecorder.fwd.hh
///
/// @brief  Forward declarations for class SilentTrajectoryRecorder.
/// @author

//  Note:   Uncomment lines beginning with "---" if you need
//          smart pointers.


#ifndef INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_fwd_hh
#define INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_fwd_hh


// External library headers
// --- #include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


// C++ headers


// Operating system headers


namespace protocols {
namespace canonical_sampling {


// Forward declarations
class SilentTrajectoryRecorder;


// Typedefs
//typedef utility::pointer::access_ptr< SilentTrajectoryRecorder > SilentTrajectoryRecorderAP;
typedef utility::pointer::shared_ptr< SilentTrajectoryRecorder > SilentTrajectoryRecorderOP;


// Smart pointer required functions
// --- void owning_ptr_acquire(SilentTrajectoryRecorder *);
// --- void owning_ptr_release(SilentTrajectoryRecorder *);

} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_SilentTrajectoryRecorder_FWD_HH
