// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/SetAACompositionPotential.fwd.hh
/// @brief  forward declaration for SetAACompositionPotential
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_potentials_SetAACompositionPotential_FWD_HH
#define INCLUDED_protocols_fldsgn_potentials_SetAACompositionPotential_FWD_HH

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace fldsgn {

/// @brief forward declaration for SetAACompositionPotential
class SetAACompositionPotential;


/// @brief SetAACompositionPotential owning pointer
typedef utility::pointer::shared_ptr< SetAACompositionPotential > SetAACompositionPotential_OP;


/// @brief SetAACompositionPotential const owning pointer
typedef utility::pointer::shared_ptr< SetAACompositionPotential const > SetAACompositionPotential_COP;


/// @brief SetAACompositionPotential access pointer
typedef utility::pointer::weak_ptr< SetAACompositionPotential > SetAACompositionPotential_AP;


/// @brief SetAACompositionPotential const access pointer
typedef utility::pointer::weak_ptr< SetAACompositionPotential const > SetAACompositionPotential_CAP;


} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_fldsgn_BluePrintBDR_FWD_HH */
