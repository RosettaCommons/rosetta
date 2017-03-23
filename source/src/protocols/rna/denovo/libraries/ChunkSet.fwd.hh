// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/libraries/ChunkSet.fwd.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_libraries_ChunkSet_FWD_HH
#define INCLUDED_protocols_farna_libraries_ChunkSet_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace libraries {

class ChunkSet;
typedef utility::pointer::shared_ptr< ChunkSet > ChunkSetOP;
typedef utility::pointer::shared_ptr< ChunkSet const > ChunkSetCOP;

} //libraries
} //denovo
} //rna
} //protocols

#endif
