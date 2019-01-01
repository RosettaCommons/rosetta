// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/rdc/RDCEnergy.fwd.hh
/// @brief   forward declaration for RDCEnergy.hh
/// @details last Modified: 08/03/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_rdc_RDCEnergy_FWD_HH
#define INCLUDED_protocols_nmr_rdc_RDCEnergy_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace nmr {
namespace rdc {

class RDCEnergy;

typedef utility::pointer::shared_ptr< RDCEnergy > RDCEnergyOP;
typedef utility::pointer::shared_ptr< RDCEnergy const > RDCEnergyCOP;
typedef utility::pointer::weak_ptr< RDCEnergy > RDCEnergyAP;
typedef utility::pointer::weak_ptr< RDCEnergy const > RDCEnergyCAP;


} // namespace rdc
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_rdc_RDCEnergy_FWD_HH
