// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/frag_picker/frag_movers/ReduceFragLengthMover
///
/// @brief      Converts a larger frag set such as 9mer to a smaller one such as 3mers
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#ifndef INCLUDED_protocols_frag_picker_frag_movers_ReduceFragLengthMover_fwd_hh
#define INCLUDED_protocols_frag_picker_frag_movers_ReduceFragLengthMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace frag_picker {
namespace frag_movers {

class ReduceFragLengthMover;
typedef utility::pointer::shared_ptr< ReduceFragLengthMover > ReduceFragLengthMoverOP;
typedef utility::pointer::shared_ptr< ReduceFragLengthMover const > ReduceFragLengthMoverCOP;

}//frag_movers
}//frag_picker
}//protocols

#endif