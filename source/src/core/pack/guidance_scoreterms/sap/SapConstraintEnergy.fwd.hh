// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintEnergy.fwd.hh
/// @brief Energy method that enforces the sap_constraint
/// @details sap_constraint constrains your protein to be soluble.
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergy_fwd_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

class SapConstraintEnergy;

typedef utility::pointer::shared_ptr< SapConstraintEnergy > SapConstraintEnergyOP;


} //sap
} //guidance_scoreterms
} //pack
} //core


#endif
