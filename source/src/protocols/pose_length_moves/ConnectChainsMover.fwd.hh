// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/pose_length_moves/ConnectChainsMover.fwd.hh
/// @brief Connects chains with native like loops.
///
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_pose_length_moves_ConnectChainsMover_fwd_hh
#define INCLUDED_protocols_pose_length_moves_ConnectChainsMover_fwd_hh

// Unit Headers

// Project headers

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace pose_length_moves {

class ConnectChainsMover;
typedef utility::pointer::shared_ptr< ConnectChainsMover > ConnectChainsMoverOP;
typedef utility::pointer::shared_ptr< ConnectChainsMover const > ConnectChainsMoverCOP;

} // pose_length_moves
} // protocols


#endif
