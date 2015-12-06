// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/farna/libraries/RNA_ChunkLibrary.hh
/// @brief  Forward declarations header
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_ChunkLibrary_fwd_hh
#define INCLUDED_protocols_rna_RNA_ChunkLibrary_fwd_hh

// ObjexxFCL Headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers

namespace protocols {
namespace farna {
namespace libraries {

class ChunkSet;
class RNA_ChunkLibrary;

typedef utility::pointer::shared_ptr< RNA_ChunkLibrary > RNA_ChunkLibraryOP;
typedef utility::pointer::shared_ptr< RNA_ChunkLibrary const > RNA_ChunkLibraryCOP;
typedef utility::pointer::shared_ptr< ChunkSet > ChunkSetOP;

} //libraries
} //farna
} //protocols

#endif
