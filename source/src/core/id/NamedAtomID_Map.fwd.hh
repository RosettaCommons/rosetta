// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/NamedAtomID_Map.fwd.hh
/// @brief  core::id::NamedAtomID_Map forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com), Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_id_NamedAtomID_Map_fwd_hh
#define INCLUDED_core_id_NamedAtomID_Map_fwd_hh


namespace core {
namespace id {


// Forward
template< typename T > class NamedAtomID_Map;

typedef  NamedAtomID_Map< bool >  NamedAtomID_Mask;

} // namespace id
} // namespace core

#endif // INCLUDED_core_id_NamedAtomID_Map_FWD_HH
