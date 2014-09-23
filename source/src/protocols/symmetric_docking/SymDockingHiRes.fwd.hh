// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/symmetric_docking/SymDockingHiRes.fwd.hh
///
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockingHiRes_fwd_hh
#define INCLUDED_protocols_symmetric_docking_SymDockingHiRes_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace symmetric_docking {


class SymDockingHiRes; // fwd declaration
typedef utility::pointer::shared_ptr< SymDockingHiRes > SymDockingHiResOP;
typedef utility::pointer::shared_ptr< SymDockingHiRes const > SymDockingHiResCOP;


} // namespace symmetric_docking
} // namespace protocols

#endif // INCLUDED_protocols_symmetric_docking_SymDockingHiRes_FWD_HH

