// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/frag_picker/fragment_movers/MergeFragSetMover
///
/// @brief      Merges two fragment sets
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#ifndef INCLUDED_protocols_frag_picker_fragment_movers_MergeFragSetMover_fwd_hh
#define INCLUDED_protocols_frag_picker_fragment_movers_MergeFragSetMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace frag_picker {
namespace frag_movers {

class MergeFragSetMover;
typedef utility::pointer::shared_ptr< MergeFragSetMover > MergeFragSetMoverOP;
typedef utility::pointer::shared_ptr< MergeFragSetMover const > MergeFragSetMoverCOP;

}//frag_movers
}//frag_picker
}//protocols

#endif
