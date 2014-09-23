// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallChunkFilter.fwd.hh
/// @brief  forward declaration for VallChunkFilter
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_VallChunkFilter_fwd_hh
#define INCLUDED_protocols_frag_picker_VallChunkFilter_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {

/// @brief forward declaration for VallChunkFilter
class VallChunkFilter;

typedef utility::pointer::shared_ptr<VallChunkFilter> VallChunkFilterOP;
typedef utility::pointer::shared_ptr<VallChunkFilter const> VallChunkFilterCOP;

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallChunkFilter_FWD_HH */
