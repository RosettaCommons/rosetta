// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/OutputMovers.fwd.hh
/// @brief  OutputMovers forward declarations header
/// @author

#ifndef INCLUDED_protocols_moves_OutputMovers_fwd_hh
#define INCLUDED_protocols_moves_OutputMovers_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

//Forwards and OP typedefs
class PDBDumpMover;
typedef utility::pointer::shared_ptr< PDBDumpMover > PDBDumpMoverOP;
typedef utility::pointer::shared_ptr< PDBDumpMover const > PDBDumpMoverCOP;

class ProfilerMover;
typedef utility::pointer::shared_ptr< ProfilerMover > ProfilerMoverOP;
typedef utility::pointer::shared_ptr< ProfilerMover const > ProfilerMoverCOP;

class MCShowMover;
typedef utility::pointer::shared_ptr< MCShowMover > MCShowMoverOP;
typedef utility::pointer::shared_ptr< MCShowMover const > MCShowMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_OutputMovers_FWD_HH
