// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.fwd.hh
/// @brief  Class declaration for the ChunkLibraryInnerLarvalJob
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_FWD_HH
#define INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace chunk_library {

class ChunkLibraryInnerLarvalJob;

typedef utility::pointer::shared_ptr< ChunkLibraryInnerLarvalJob > ChunkLibraryInnerLarvalJobOP;
typedef utility::pointer::shared_ptr< ChunkLibraryInnerLarvalJob const > ChunkLibraryInnerLarvalJobCOP;

typedef utility::vector1< ChunkLibraryInnerLarvalJobOP > ChunkLibraryInnerLarvalJobs;

} // namesapce chunk_library
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_FWD_HH
