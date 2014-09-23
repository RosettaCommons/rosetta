// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/SetSecStructEnergies.fwd.hh
/// @brief  forward declaration for SetSecStructEnergies
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_potentials_SetSecStructEnergies_FWD_hh
#define INCLUDED_protocols_fldsgn_potentials_SetSecStructEnergies_FWD_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {

/// @brief forward declaration for SetSecStructEnergies
class SetSecStructEnergies;


/// @brief SetSecStructEnergies owning pointer
typedef utility::pointer::shared_ptr< SetSecStructEnergies > SetSecStructEnergies_OP;


/// @brief SetSecStructEnergies const owning pointer
typedef utility::pointer::shared_ptr< SetSecStructEnergies const > SetSecStructEnergies_COP;


/// @brief SetSecStructEnergies access pointer
typedef utility::pointer::weak_ptr< SetSecStructEnergies > SetSecStructEnergies_AP;


/// @brief SetSecStructEnergies const access pointer
typedef utility::pointer::weak_ptr< SetSecStructEnergies const > SetSecStructEnergies_CAP;

} // namespace potentials
} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_fldsgn_BluePrintBDR_FWD_HH */
