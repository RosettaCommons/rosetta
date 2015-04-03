// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/DunbrackRotamer.fwd.hh
/// @brief  Typedefs and forward declarations for class DunbrackRotamer
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_ChiSet_fwd_hh
#define INCLUDED_core_pack_dunbrack_ChiSet_fwd_hh

// Package headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace pack {
namespace dunbrack {

class ChiSet;
typedef utility::pointer::shared_ptr< ChiSet > ChiSetOP;
typedef utility::pointer::shared_ptr< ChiSet const > ChiSetCOP;

}
}
}

#endif //INCLUDED_core_pack_dunbrack_ChiSet_FWD_HH


