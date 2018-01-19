// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.fwd.hh
/// @brief A helper class to assist in parsing rotamer libraries.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_fwd_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace dunbrack {

class RotamericSingleResidueDunbrackLibraryParser;

typedef utility::pointer::shared_ptr< RotamericSingleResidueDunbrackLibraryParser > RotamericSingleResidueDunbrackLibraryParserOP;
typedef utility::pointer::shared_ptr< RotamericSingleResidueDunbrackLibraryParser const > RotamericSingleResidueDunbrackLibraryParserCOP;

} //core
} //pack
} //dunbrack

#endif //INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_fwd_hh
