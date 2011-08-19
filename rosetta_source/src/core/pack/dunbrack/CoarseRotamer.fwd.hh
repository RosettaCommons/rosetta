// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/CoarseSingleResidueLibrary.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_CoarseRotamer_fwd_hh
#define INCLUDED_core_pack_dunbrack_CoarseRotamer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {

class CoarseRotamer;
typedef utility::pointer::owning_ptr< CoarseRotamer > CoarseRotamerOP;
typedef utility::pointer::owning_ptr< CoarseRotamer const > CoarseRotamerCOP;

class CoarseRotamerSet;
typedef utility::pointer::owning_ptr< CoarseRotamerSet > CoarseRotamerSetOP;
typedef utility::pointer::owning_ptr< CoarseRotamerSet const > CoarseRotamerSetCOP;


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_CoarseSingleResidueLibrary_FWD_HH
