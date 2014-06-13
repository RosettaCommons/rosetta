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


#ifndef INCLUDED_core_pack_dunbrack_DunbrackRotamer_fwd_hh
#define INCLUDED_core_pack_dunbrack_DunbrackRotamer_fwd_hh

// Package headers
// AUTO-REMOVED #include <core/scoring/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/types.hh>
#include <utility/vector1.fwd.hh>



namespace core {
namespace pack {
namespace dunbrack {

/// Low precision in Dunbrack rotamer library suggests no need to store
/// chi dihedrals and standard deviations with 64 bits.
typedef float DunbrackReal;

typedef utility::vector1< Size > RotVector;
typedef utility::vector1< Real > ChiVector;
typedef utility::vector1< Real > AngleVector;

Size const ONE = 1;
Size const TWO = 2;
Size const THREE = 3;
Size const FOUR = 4;

/// @brief forward declaration; default precision is DunbrackReal
template < Size S, class P = DunbrackReal >
class DunbrackRotamerMeanSD;

/// @brief forward declaration; default precision is DunbrackReal;
template < Size S, class P = DunbrackReal >
class PackedDunbrackRotamer;

/// @brief forward declaration; default precision is DunbrackReal
template < Size S, class P = DunbrackReal >
class DunbrackRotamer;

class RotamerBuildingData;

class DunbrackRotamerSampleData;

/// DOUG DOUG DOUG
template < Size T >
class RotamericData;

} // dunbrack
} // pack
} // core

#endif
