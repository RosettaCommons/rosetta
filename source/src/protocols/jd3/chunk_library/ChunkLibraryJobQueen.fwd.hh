// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/ChunkLibraryJobQueen.fwd.hh
/// @brief  class declaration for ChunkLibraryJobQueen
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_chunk_library_ChunkLibraryJobQueen_FWD_HH
#define INCLUDED_protocols_jd3_chunk_library_ChunkLibraryJobQueen_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace chunk_library {

class ChunkLibraryJobQueen;
class PreliminaryLarvalJob;

typedef utility::pointer::shared_ptr< ChunkLibraryJobQueen > ChunkLibraryJobQueenOP;
typedef utility::pointer::shared_ptr< ChunkLibraryJobQueen const > ChunkLibraryJobQueenCOP;

typedef utility::pointer::shared_ptr< PreliminaryLarvalJob > PreliminaryLarvalJobOP;
typedef utility::pointer::shared_ptr< PreliminaryLarvalJob const > PreliminaryLarvalJobCOP;

} // namespace chunk_library
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_chunk_library_ChunkLibraryJobQueen_FWD_HH
