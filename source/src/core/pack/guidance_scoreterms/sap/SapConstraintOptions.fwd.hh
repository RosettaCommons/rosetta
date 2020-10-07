// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh
/// @brief Options for SapConstraint
/// @details Contains all the options and a method to transform sap_score into sap_constraint
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintOptions_fwd_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintOptions_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

class SapConstraintOptions;

typedef utility::pointer::shared_ptr< SapConstraintOptions > SapConstraintOptionsOP;
typedef utility::pointer::shared_ptr< SapConstraintOptions const > SapConstraintOptionsCOP;

} //sap
} //guidance_scoreterms
} //pack
} //core


#endif
