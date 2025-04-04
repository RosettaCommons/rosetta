// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMBondAngleLibrary.fwd.hh
/// @brief  core::scoring::mm::MMBondAngleLibrary forward declarations
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleLibrary_fwd_hh
#define INCLUDED_core_scoring_mm_MMBondAngleLibrary_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace mm {


// Forward
class MMBondAngleLibrary;

typedef  utility::pointer::weak_ptr< MMBondAngleLibrary > MMBondAngleLibraryAP;
typedef  utility::pointer::weak_ptr< MMBondAngleLibrary const > MMBondAngleLibraryCAP;
typedef  utility::pointer::shared_ptr< MMBondAngleLibrary > MMBondAngleLibraryOP;
typedef  utility::pointer::shared_ptr< MMBondAngleLibrary const > MMBondAngleLibraryCOP;

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_bondangle_library_FWD_HH
