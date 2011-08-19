// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/Chunk.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_CHUNK_FWD_HH_
#define PROTOCOLS_NONLOCAL_CHUNK_FWD_HH_

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace nonlocal {

class Chunk;
typedef utility::pointer::owning_ptr<Chunk> ChunkOP;
typedef utility::pointer::owning_ptr<Chunk const> ChunkCOP;

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_CHUNK_FWD_HH_
