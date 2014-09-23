// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallChunk.fwd.hh
/// @brief  forward declaration for VallChunk
/// @author Dominik Gront (dgront@gmail.com)

#ifndef INCLUDED_protocols_frag_picker_VallChunk_fwd_hh
#define INCLUDED_protocols_frag_picker_VallChunk_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {

/// @brief forward declaration for VallChunk
class VallChunk;

typedef utility::pointer::owning_ptr<VallChunk> VallChunkOP;
typedef utility::pointer::owning_ptr<VallChunk const> VallChunkCOP;
typedef utility::pointer::access_ptr<VallChunk> VallChunkAP;
typedef utility::pointer::access_ptr<VallChunk const> VallChunkCAP;

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallChunk_FWD_HH */
