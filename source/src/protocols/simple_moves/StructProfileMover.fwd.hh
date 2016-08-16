// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/StructProfileMover.fwd.hh
/// @brief Quickly generates a structure profile.
///
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_simple_moves_StructProfileMover_fwd_hh
#define INCLUDED_protocols_simple_moves_StructProfileMover_fwd_hh

// Unit Headers
// you cannot #include yourself #include <protocols/moves/StructProfileMover.fwd.hh>

// Project headers

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_moves {

class StructProfileMover;
typedef utility::pointer::shared_ptr< StructProfileMover > StructProfileMoverOP;
typedef utility::pointer::shared_ptr< StructProfileMover const > StructProfileMoverCOP;

} // simple_moves
} // protocols


#endif

