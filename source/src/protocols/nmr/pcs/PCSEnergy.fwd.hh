// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSEnergy.fwd.hh
/// @brief   forward declaration for PCSEnergy.hh
/// @details last Modified: 07/19/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pcs_PCSEnergy_FWD_HH
#define INCLUDED_protocols_nmr_pcs_PCSEnergy_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace nmr {
namespace pcs {

class PCSEnergy;

typedef utility::pointer::shared_ptr< PCSEnergy > PCSEnergyOP;
typedef utility::pointer::shared_ptr< PCSEnergy const > PCSEnergyCOP;
typedef utility::pointer::weak_ptr< PCSEnergy > PCSEnergyAP;
typedef utility::pointer::weak_ptr< PCSEnergy const > PCSEnergyCAP;


} // namespace pcs
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pcs_PCSEnergy_FWD_HH
