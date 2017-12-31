// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library_inputters/PDBChunkLibraryInputter.fwd.hh
/// @brief  Forward declaration of the %PDBChunkLibraryInputter class for initializing Poses from .pdb or .pdb.gz files
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_jd3_chunk_library_inputters_PDBChunkLibraryInputter_FWD_HH
#define INCLUDED_protocols_jd3_chunk_library_inputters_PDBChunkLibraryInputter_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace chunk_library_inputters {

class PDBChunkLibraryInputter;

typedef utility::pointer::shared_ptr< PDBChunkLibraryInputter > PDBChunkLibraryInputterOP;
typedef utility::pointer::shared_ptr< PDBChunkLibraryInputter const > PDBChunkLibraryInputterCOP;


} // namespace chunk_library_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_PDBChunkLibraryInputter_FWD_HH
