// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/frag_picker/fragment_movers/FragSetFromH5Mover
///
/// @brief      Pick fragments from the H5 database.
/// @details    Allows based on ABEGO def,SSdef,PredictedSS,9mers,3mers pushes to stack
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#ifndef INCLUDED_protocols_frag_picker_frag_movers_FragSetFromH5Mover_fwd_hh
#define INCLUDED_protocols_frag_picker_frag_movers_FragSetFromH5Mover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

#include <protocols/frag_picker/frag_movers/FragSetFromH5Mover.hh>

namespace protocols {
namespace frag_picker {
namespace frag_movers {

class FragSetFromH5Mover;
typedef utility::pointer::shared_ptr< FragSetFromH5Mover > FragSetFromH5MoverOP;
typedef utility::pointer::shared_ptr< FragSetFromH5Mover const > FragSetFromH5MoverCOP;

}//frag_movers
}//frag_picker
}//protocols

#endif