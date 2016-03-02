// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/NearNativeLoopCloser.fwd.hh
/// @brief Uses KIC/CCD and a fast Look back to close loops near native
///
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_pose_length_moves_NearNativeLoopCloser_fwd_hh
#define INCLUDED_protocols_pose_length_moves_NearNativeLoopCloser_fwd_hh

// Unit Headers
// you cannot #include yourself #include <protocols/moves/NearNativeLoopCloser.fwd.hh>

// Project headers

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace pose_length_moves {

class PossibleLoop;
typedef utility::pointer::shared_ptr< PossibleLoop > PossibleLoopOP;
typedef utility::pointer::shared_ptr< PossibleLoop const > PossibleLoopCOP;

class NearNativeLoopCloser;
typedef utility::pointer::shared_ptr< NearNativeLoopCloser > NearNativeLoopCloserOP;
typedef utility::pointer::shared_ptr< NearNativeLoopCloser const > NearNativeLoopCloserCOP;

} // pose_length_moves
} // protocols


#endif
