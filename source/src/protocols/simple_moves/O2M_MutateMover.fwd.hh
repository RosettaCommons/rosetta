// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/O2M_MutateMover.fwd.hh
/// @brief  O2M_MutateMover forward declarations header
/// @author Ken Jung

#ifndef INCLUDED_protocols_simple_moves_O2M_MutateMover_fwd_hh
#define INCLUDED_protocols_simple_moves_O2M_MutateMover_fwd_hh

// Utility headers
#include <boost/shared_ptr.hpp>

namespace protocols {
namespace simple_moves {

class O2M_MutateMover;
typedef boost::shared_ptr<O2M_MutateMover> O2M_MutateMoverSP;

}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_O2M_MutateMover_fwd_hh
