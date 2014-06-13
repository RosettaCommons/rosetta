// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMTorsionLibrary.fwd.hh
/// @brief  core::scoring::mm::MMTorsionLibrary forward declarations
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMTorsionLibrary_fwd_hh
#define INCLUDED_core_scoring_mm_MMTorsionLibrary_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace scoring {
namespace mm {


// Forward
class MMTorsionLibrary;

typedef  utility::pointer::access_ptr< MMTorsionLibrary > MMTorsionLibraryAP;
typedef  utility::pointer::access_ptr< MMTorsionLibrary const > MMTorsionLibraryCAP;
typedef  utility::pointer::owning_ptr< MMTorsionLibrary > MMTorsionLibraryOP;
typedef  utility::pointer::owning_ptr< MMTorsionLibrary const > MMTorsionLibraryCOP;

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_torsion_library_FWD_HH
