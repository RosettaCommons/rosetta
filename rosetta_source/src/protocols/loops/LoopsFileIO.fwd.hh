// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopsFileIO.fwd.hh
/// @brief  LoopsFileIO forward declarations header
/// @author Brian D. Weitzner


#ifndef INCLUDED_protocols_loops_LoopsFileIO_FWD_HH
#define INCLUDED_protocols_loops_LoopsFileIO_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {

// Forward
class LoopsFileIO;

typedef utility::pointer::owning_ptr< LoopsFileIO > LoopsFileIOOP;
typedef utility::pointer::owning_ptr< LoopsFileIO const > LoopsFileIOCOP;

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_LoopsFileIO_FWD_HH
